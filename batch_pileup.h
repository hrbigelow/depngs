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
    khint_t indel_key;
    unsigned ct;
};


void batch_pileup_init(unsigned n_samples);
void batch_pileup_free();

void reference_seq_init(char *fasta_file);
void reference_seq_free();


/* perform entire tally phase, for basecalls and for indels for one
   sample. update tls.tally_end to reflect furthest position seen. */
void tally_pileup_stats(const struct managed_buf bam, unsigned s);

/* summarize statistics for sample s, up to tls.tally_end, or to
   completion if there is no more input available */
void summarize_pileup_stats(unsigned s);


/* provide basecall stats for a sample at current position, or the
   null statistic. */
struct base_count pileup_basecall_stats(unsigned s);

/* provide (b,q,s) stats for a sample at current position.  *cts is
   reallocated as necessary. *n_cts set to number of distinct stats
   that are populated. */
void pileup_bqs_stats(unsigned s, struct bqs_count **cts, unsigned *n_cts);

/* provide indel stats for a sample at current position */
void pileup_indel_stats(unsigned s, struct indel_count **cts, unsigned *n_cts);

/* produce pileup strings from current position for sample s,
   storing in call and qual, respectively */
void pileup_strings(unsigned s, struct managed_buf *call, struct managed_buf *qual);

/* advance internal iterator position for base counts */
int next_base();

/* advance internal iterator position for indel counts */
int next_indel();

/* return name of reference at current position */
char *pileup_current_refname();

/* return reference base at current position */
char pileup_current_refbase();

/* return the current base position */
unsigned pileup_current_pos();

/* clear statistics that are no longer needed. call this function
   after next_base() and next_indel() return 1. */
void clear_finished_stats();


#endif /* _BATCH_PILEUP_H */
