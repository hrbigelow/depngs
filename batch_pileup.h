#ifndef _BATCH_PILEUP_H
#define _BATCH_PILEUP_H

/*

  Initialization:
  
  1. Call batch_pileup_init() and reference_seq_init().
  
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
#include "khash.h"
#include "cache.h"

struct base_count {
    unsigned ct[4];
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


/* called once for each thread when starting to use the batch_pileup
   functionality. */
void
batch_pileup_thread_init(unsigned n_samples);

/* called once for each thread at end of using batch_pileup
   functionality */
void batch_pileup_thread_free();

/* called once at start of program */
void reference_seq_init(char *fasta_file);
void reference_seq_free();


/* perform entire tally phase, for basecalls and for indels for one
   sample. update tls.tally_end to reflect furthest position seen. */
void tally_pileup_stats(const struct managed_buf bam, unsigned s);

/* summarize statistics for sample s, up to tls.tally_end, or to
   completion if there is no more input available.  initializes
   base_ct, base_cur, and base_end, as well as indel_ct, indel_cur,
   and indel_end.  base_cur and indel_cur are set to the beginning. */
void summarize_pileup_stats(unsigned s);


/* provide basecall stats for a sample at current position, or the
   null statistic. */
struct base_count
pileup_basecall_stats(unsigned s);

/* provide (b,q,s) stats for a sample at current position.  *cts is
   reallocated as necessary. *n_cts set to number of distinct stats
   that are populated. */
void
pileup_bqs_stats(unsigned s, struct bqs_count **cts, unsigned *n_cts);

/* provide indel stats for a sample at current position */
void
pileup_indel_stats(unsigned s, struct indel_count **cts, unsigned *n_cts);

/* produce pileup strings from current position for sample s,
   storing in call and qual, respectively */
void pileup_strings(unsigned s, struct managed_buf *call, struct managed_buf *qual);


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


struct pileup_data {
    struct managed_buf calls, quals;
    unsigned read_depth; /* total read depth */
    unsigned used_read_depth; /* read depth used (due to above-quality
                                 threshold) */
};


/* produce pileup data (calls, quals, and read depths) from current
   position for sample s */
void
pileup_current_data(unsigned s, struct pileup_data *pd);


/* clear statistics that are no longer needed. call this function
   after next_base() and next_indel() return 1. */
void pileup_clear_finished_stats();


#endif /* _BATCH_PILEUP_H */
