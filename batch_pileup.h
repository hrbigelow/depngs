#ifndef _BATCH_PILEUP_H
#define _BATCH_PILEUP_H

/*
  Initialization:
  
  1. Call batch_pileup_thread_init() once for each thread.  Call
     batch_pileup_init() once for the entire program.
  
  Use:
  
  1. Obtain an array of blocks of BAM records that overlap a given
  region, then call batch_tally() to compile all statistics.
  
  2. Call basecall_stats(), bqs_stats(), and indel_stats() to get
  views of these summary statistics for the current base or indel
  iterator position.
  
  3. Call next_base() and next_indel() to advance base and indel
  iterators and repeat from step 2. These return 1 if end is
  reached.
  
  4. Once you reach the end, call clear_finished_stats() to release
  the set of stats that are no longer needed.  (Residual
  statistics that are only partially compiled will be retained,
  and added to in the next call to batch_tally().
  
  Cleanup:
  
  6. Call batch_pileup_free() and reference_seq_free();
  
 */
#include <stdint.h>
#include "htslib/sam.h"
#include "khash.h"
#include "cache.h"
#include "compat_util.h"

struct indel_seq {
    char is_ins;
    char seq[FLEX_ARRAY]; /* zero-terminated */
};


struct base_count {
    unsigned ct_filt[4]; /* # CIGAR match reads by base call, with
                              qual >= min_quality_score */
    unsigned n_match_lo_q; /* # CIGAR match reads with qual < min_quality_score */
    unsigned n_match_hi_q; /* # CIGAR match reads with qual >= min_quality_score */
};


struct bqs_count {
    uint32_t base: 2;
    uint32_t qual: 7;
    uint32_t strand: 1;
    uint32_t ct: 22; /* 4,194,304 */
};

struct indel_count {
    khiter_t indel_itr;
    unsigned ct;
};


struct indel_pair_count {
    unsigned indel_id; /* can be used to find the indel */
    unsigned count[2];
};

/* called once for each thread when starting to use the batch_pileup
   functionality. */
void
batch_pileup_thread_init(unsigned n_samples);

/* called once for each thread at end of using batch_pileup
   functionality */
void
batch_pileup_thread_free();

/* initialize all program-wide static data (refseq and
   min_quality_score) */
void
batch_pileup_init(const char *fasta_file, unsigned min_qual);

void
batch_pileup_free();


/* perform entire tally phase, for basecalls and for indels for one
   sample. update tls.tally_end to reflect furthest position seen. */
void
tally_pileup_stats(const struct managed_buf bam, unsigned s);

/* summarize statistics for sample s, up to tls.tally_end, or to
   completion if there is no more input available.  initializes
   base_ct, base_cur, and base_end, as well as indel_ct, indel_cur,
   and indel_end.  base_cur and indel_cur are set to the beginning. */
void
summarize_pileup_stats(unsigned s);


/* return an initialized base_count structure for sample s at the
   current position.  see 'struct base_count' for details. */
struct base_count
pileup_basecall_stats(unsigned s);


/* provide (b,q,s) stats for a sample at current position.  *cts is
   reallocated as necessary. *n_cts set to number of distinct stats
   that are populated. */
void
pileup_bqs_stats(unsigned s, struct bqs_count **cts, unsigned *n_cts);

/* provide indel stats for a sample at current position.  They will be
   sorted ascending by ict.indel_itr field */
void
pileup_indel_stats(unsigned s, struct indel_count **cts, unsigned *n_cts);

/* provide pairwise indel stats.  reallocates *pair_cts as
   necessary */
void
pileup_pair_indel_stats(struct indel_count *cts1, unsigned n_cts1,
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
int pileup_next_pos();

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
free_pileup_data(struct pileup_data *pd);

/* produce pileup data (calls, quals, and read depths) from current
   position for sample s */
void
pileup_current_data(unsigned s, struct pileup_data *pd);


/* clear statistics that are no longer needed. call this function
   after next_base() and next_indel() return 1. */
void
pileup_clear_finished_stats();


#endif /* _BATCH_PILEUP_H */
