#ifndef _CHUNK_STRATEGY_H
#define _CHUNK_STRATEGY_H

#include "locus_range.h"

/* tallys statistics that allow the reader to adjust its strategy
   during program execution. one of two strategies is used to estimate
   the number of bytes left.

   1.  if do_range_estimation is set, we have been given a list of
       locus ranges to process.  n_loci_total is initialized at the
       start. n_loci_read, n_bytes_read are maintained.  A running
       estimate of bytes per locus is maintained, and from this, an
       estimate of bytes left.

   2.  otherwise, we are processing the whole file.  n_bytes_total and
       n_bytes_read are used to directly account for the number of
       bytes left.

       
 */
struct chunk_strategy {
    /* marker that informs all threads where to resume reading */
    struct contig_pos pos;

    unsigned n_files;
    unsigned long *n_bytes_read;
    unsigned long n_loci_total;
    unsigned long n_loci_read;

    unsigned long max_bytes_small_chunk; /* # bytes left to switch to
                                              small chunk strategy */
    unsigned long small_chunk_size;
    unsigned long default_bytes_per_locus;
};


/* call once at start of program */
void
chunk_strategy_init(unsigned n_files);


/* call once at end of program */
void
chunk_strategy_free();


/* call if you are re-running a new chunk of input that is within a
  thread_queue_run() call. */
void
chunk_strategy_reset(unsigned long n_loci);



/* call this if a range file is given */
void cs_init_by_range(unsigned n_loci_total, unsigned n_files);

/* call this if no range file is given */
void cs_init_whole_file(unsigned n_files);

/* call to reset the position, if you want to re-process input. */
void
cs_stats_reset_pos();


/* configure the chunking strategy.  if < max_bytes_small_chunk of
   input remain, switch to small chunks of size small_chunk_size.  if
   we are doing range estimation, use default_bytes_per_locus as an
   initial estimate (before any input is read) in order to convert
   from a locus count estimate to a bytes estimate.
 */
void cs_set_defaults(unsigned long max_bytes_small_chunk,
                     unsigned long small_chunk_size,
                     unsigned long default_bytes_per_locus);


/* estimate the bytes wanted based on the strategy */
unsigned long
cs_get_bytes_wanted(unsigned n_files); 


#endif /* _CHUNK_STRATEGY_H */
