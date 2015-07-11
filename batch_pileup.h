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

void batch_pileup_init(unsigned n_samples);
void batch_pileup_free();

void reference_seq_init(char *fasta_file);
void reference_seq_free();


/* perform entire tally phase, for basecalls and for indels. bam[s] is
   block of bam records for sample s */
void batch_tally(const struct managed_buf *bam);

/* provide basecall stats for a sample at current position, or the
   null statistic. */
struct base_count basecall_stats(unsigned s);

/* provide (b,q,s) stats for a sample at current position.  *cts is
   reallocated as necessary. *n_cts set to number of distinct stats
   that are populated. */
void bqs_stats(unsigned s, struct bqs_count **cts, unsigned *n_cts);

/* provide indel stats for a sample at current position */
void indel_stats(unsigned s, struct indel_count **cts, unsigned *n_cts);

/* advance internal iterator position for base counts */
int next_base();

/* advance internal iterator position for indel counts */
int next_indel();

/* clear statistics that are no longer needed. call this function
   after next_base() and next_indel() return 1. */
void clear_finished_stats();


#endif /* _BATCH_PILEUP_H */
