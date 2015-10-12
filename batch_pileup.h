#ifndef _BATCH_PILEUP_H
#define _BATCH_PILEUP_H

/*
  Synopsis:

  1. At program start: batch_pileup_init()
  2. At each thread start: batch_pileup_thread_init()
  3. After reading a chunk of data:
     a. call pileup_load_refseq_ranges() once.
     b. call pileup_tally_stats() once for each sample.
     c. if there is no more input after this, call pileup_final_input()
     d. (optional) call pileup_prepare_basecalls() once for each sample
     e. (optional) call pileup_prepare_bqs() once for each sample
     f. (optional) call pileup_prepare_indels() once for each sample
  4. Main loop:
     a. call pileup_next_pos()
     b. (optional) call pileup_current_basecalls()
     c. (optional) call pileup_current_bqs()
     d. (optional) call pileup_current_indels()

  5. After main loop:
     a. call pileup_clear_stats()

  6. At end of input:
     a. At each thread end: batch_pileup_thread_free()
     b. At program end: batch_pileup_free()

  7. (optional) call pileup_reset_pos() to reset the position marker.
  This enables re-running the main loop (step 4) on the same input.

  Notes:

  The optional pileup_prepare_* functions must be called for each
  matching pileup_current_* function that you want to call.

  pileup_next_pos() advances the current position marker.  The current
  position marker (tls.cur_pos) determines the output of each
  pileup_current_* function.
  
*/

#include <stdint.h>
#include "htslib/sam.h"
#include "khash.h"
#include "cache.h"
#include "compat_util.h"
#include "bam_reader.h"


 /* n_match_lo_q + n_match_hi_q + n_match_fuzzy equal the total number
    of CIGAR match bases at this locus. ct_filt contain the
    n_match_hi_q calls broken down by basecall. (A 'fuzzy' call is a
    basecall that is not one of A,C,G, or T).  */
 struct base_count {
     unsigned ct_filt[4]; /* # CIGAR match reads by base call, with
                             qual >= min_quality_score */
     unsigned n_match_lo_q; /* # CIGAR match reads with qual < min_quality_score */
     unsigned n_match_hi_q; /* # CIGAR match reads with qual >= min_quality_score */
     unsigned n_match_fuzzy; /* # CIGAR match reads with ambiguous base call.   */
 };


struct bqs_count {
    uint32_t base: 4; /* 4-bit BAM encoding of base.  use hts.h:
                         seq_nt16_str[base] for the letter-equivalent. */
    uint32_t qual: 7;
    uint32_t strand: 1;
    uint32_t ct: 20; /* 1,048,576 */
};


struct indel {
    char is_ins;
    unsigned length;
    const char *seq;
};

struct indel_count {
    struct indel idl;
    unsigned ct;
};

/* the reified version of struct indel. seq is populated either from
   the str_hash (if it is an insertion) or by querying the reference
   genome (if deletion) */
struct indel_seq {
    char is_ins;
    char seq[FLEX_ARRAY]; /* zero-terminated */
};


/* value for the hash */
struct indel_count_ary {
    struct indel_count *i;
    unsigned n, m;
};


struct indel_pair_count {
    struct indel indel;
    unsigned count[2];
};

/* called once for each thread when starting to use the batch_pileup
   functionality. */
void
batch_pileup_thread_init(unsigned n_samples, const char *fasta_file);

/* called once for each thread at end of using batch_pileup
   functionality */
void
batch_pileup_thread_free();


struct batch_pileup_params {
    unsigned skip_empty_loci; /* if 1, pileup_next_pos() advances past
                                 loci that have no data. */
    unsigned pseudo_depth;
    unsigned min_clash_qual; /* for doubly sequenced basecalls that
                                clash to be both set to zero qual,
                                both must be above this threshold. */
};

/* skip_empty_loci: if 0, pileup_next_pos() will visit all loci in the
   region of interest.  if set to 1, skip those loci in the region of
   interest that have no data for any sample.  setting this to 0
   allows the client to account for these no-data loci in the main
   loop. */
void
batch_pileup_init(struct bam_filter_params _bam_filter,
                  struct batch_pileup_params bp_par);


void
batch_pileup_free();



/* load specific ranges of reference sequence. [qbeg, qend) defines
   the total set of (non-overlapping) ranges to consider.  subset
   defines the overlapping intersection of these ranges that will be
   loaded into tls.refseqs.  assume that each interval in [qbeg, qend)
   is on one contig, but that 'subset' may span multiple contigs. */
void
pileup_load_refseq_ranges(struct bam_scanner_info *bsi);



/* perform entire tally phase, for basecalls and for indels for one
   sample. update tls.tally_end to reflect furthest position seen. */
void
pileup_tally_stats(const struct managed_buf bam, 
                   struct bam_scanner_info *bsi,
                   unsigned s);


/* prepares internal data structures for fast retrieval.  must call
   this function before calling pileup_current_basecalls */
void
pileup_prepare_basecalls();

/* return an initialized base_count structure for sample s at the
   current position.  see 'struct base_count' for details. */
struct base_count
pileup_current_basecalls(unsigned s);


/* prepare internal data structures for fast retrieval of bqs
   information by position.  must call this function before calling
   pileup_current_bqs. */
void
pileup_prepare_bqs(unsigned s);

 /* provide (b,q,s) stats for a sample at current position.  *cts is
    reallocated as necessary. *n_cts set to number of distinct stats
    that are populated. */
void
pileup_current_bqs(unsigned s, struct bqs_count **cts, unsigned *n_cts);


/* prepare internal data structures for fast retrieval of indel
   information.  must call this function before calling
   pileup_current_indels. */
void
pileup_prepare_indels(unsigned s);

/* provide indel stats for a sample at current position.  The order of
   elements in (*cts) is defined by pos_iter_indel_count_less.  For
   indels at a current position, the order is: deletions first, then
   insertions, then by ascending length of deletion (or ascending by
   insertion iterator key). Admittedly, a weird ordering... */
void
pileup_current_indels(unsigned s, struct indel_count **cts, unsigned *n_cts);

/* merge information in indel counts 1 and 2, producing counts of
   pairs based on indel type. */
void
pileup_make_indel_pairs(struct indel_count *cts1, unsigned n_cts1,
                        struct indel_count *cts2, unsigned n_cts2,
                        struct indel_pair_count **pair_cts, unsigned *n_pair_cts);


/* populate and possibly re-allocate indel with a copy of the indel
   identified by indel_id.  caller must free *indel */
void
pileup_get_indel(unsigned indel_id, struct indel_seq **indel);


/* produce pileup strings from current position for sample s,
   storing in call and qual, respectively */
void
pileup_strings(unsigned s, struct managed_buf *call, struct managed_buf *qual);


/* advance tls.cur_pos to next position for which at least one sample
   has base or indel entries. update base_cur and indel_cur pointers
   for each sample so that calls to pileup_{basecall,bqs,indel}_stats
   return the right values. return 1 if there is a next position, 0 if
   reached the end. */
int
pileup_next_pos();

#define REFNAME_MAXLEN 300

struct pileup_locus_info {
    char refname[REFNAME_MAXLEN + 1];
    char refbase;
    unsigned pos;
};


/* returns information about current locus (not pertaining to any
   specific sample) */
void
pileup_current_info(struct pileup_locus_info *pli);


/* like pileup_current_info but avoids querying for the contig name */
struct contig_pos
pileup_current_pos();


/* description of reads data at a single-base locus position.  for non
   sample data-dependent information, see pileup_locus_info. */
struct pileup_data {
    struct managed_buf calls, quals;
    unsigned n_match_lo_q; /* # reads with CIGAR match state and qual
                                   < min_quality_score */
    unsigned n_match_hi_q; /* # reads with CIGAR match state and qual
                                   >= min_quality_score */
    unsigned n_indel; /* # reads with CIGAR indel state */
};


void
init_pileup_data(struct pileup_data *pd);


void
free_pileup_data(struct pileup_data *pd);


/* reify the indel, allocating and returning an indel_seq, using the
   current position.  Note: this function doesn't actually need the
   current position in order to retrieve sequences for insertions, but
   it does for deletions. */
struct indel_seq *
pileup_current_indel_seq(struct indel *idl);


/* produce pileup data (calls, quals, and read depths) from current
   position for sample s.  must call pileup_prepare_bqs and
   pileup_prepare_indels before calling this. */
void
pileup_current_data(unsigned s, struct pileup_data *pd);


/* clear all statistics in preparation for a new chunk */
void
pileup_clear_stats();


/* call to signal that there will be no more input.  After calling
   this, a call to summarize_pileup_stats() will summarize all
   remaining residual statistics. */
void
pileup_final_input();


/* call to reset the position marker to its beginning state.  This
   allows you to reprocess the same input in a main loop using
   pileup_next_pos() */
void
pileup_reset_pos();

/* retrieve the bam record filtering parameters in effect */
struct bam_filter_params
pileup_get_filter_params();

#endif /* _BATCH_PILEUP_H */
