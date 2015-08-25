#ifndef _CHUNK_STRATEGY_H
#define _CHUNK_STRATEGY_H

#include "locus_range.h"

/* The chunking strategy attempts to prevent thread starvation as the
   program nears the end of the input.  The input for a given sample
   of total size T, is divided into three zones of decreasing size.
   The first zone is of size T - S2 - S1.  The second and third are of
   sizes S2 and S1.  S1 and S2 are given by the user, where S1 < S2.
   T, however, is estimated for each sample, based on the known total
   number of loci, and the estimated number of bytes per locus.

   If T is less than (S1 + S2), then the first zone doesn't exist.  In
   each zone, the function cs_max_bytes_wanted returns the size of the
   zone divided by n_threads.  This ensures that at least n_threads
   will get a piece of that zone.

   The difficulty in this is that we cannot simply request input based
   on a total number of loci wanted, because this violates a hard
   memory constraint. */
struct chunk_strategy {

    /* global configuration defining input and strategy */
    unsigned n_files, n_threads;
    const struct contig_region *query_regions;
    unsigned n_query_regions;
    unsigned long bytes_zone2, bytes_zone3;
 
    /* configuration that can be refreshed when new span is set */
    struct contig_span total_span;
    unsigned long n_total_loci;
    
    /* running per-sample statistics needed for estimating bytes
       left. */
    unsigned long *n_all_bytes_read;
    unsigned long n_all_loci_read;

    /* current position for starting the next read */
    struct contig_pos cur_pos;
};

extern struct chunk_strategy cs_stats;


/* call once at start of program. */
void
chunk_strategy_init(unsigned n_files, unsigned n_threads,
                    const char *locus_range_file,
                    const char *fasta_file,
                    unsigned long bytes_zone2,
                    unsigned long bytes_zone3);


/* call once at end of program */
void
chunk_strategy_free();


/* call if you are re-running a new chunk of input that is within a
  thread_queue_run() call. resets everything except n_all_loci_read,
  n_all_bytes_read. */
void
chunk_strategy_reset();


/* call to prepare for processing a new span. */
void
chunk_strategy_set_span(struct contig_span span);


/* estimate the bytes wanted based on the strategy */
unsigned long
cs_max_bytes_wanted(); 


#endif /* _CHUNK_STRATEGY_H */
