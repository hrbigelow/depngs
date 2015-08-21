#ifndef _CHUNK_STRATEGY_H
#define _CHUNK_STRATEGY_H

#include "locus_range.h"

/* Provides a global resource (to be accessed by one thread at a time)
   to guide which ranges of input to process.  query_regions define
   the subset of genomic loci to consider, and remains constant
   through the life of the program.  span defines the current subset
   of these ranges to load.  The client will update span as it reads
   input.  Other fields provide a mechanism for dividing up the
   remaining input to keep all threads occupied.
 */
struct chunk_strategy {

    unsigned n_files;
    unsigned n_threads;
    unsigned long *n_all_bytes_read;
    unsigned long n_all_loci_read;
    unsigned long n_loci_total;
    unsigned long n_loci_read;
    const struct contig_region *query_regions;
    unsigned n_query_regions;
    unsigned long n_min_absolute_bytes; /* min for cs_max_bytes_wanted() */
    struct contig_span span;
};

extern struct chunk_strategy cs_stats;


/* call once at start of program */
void
chunk_strategy_init(unsigned n_files, unsigned n_threads,
                    unsigned long n_min_absolute_bytes,
                    const char *locus_range_file,
                    const char *fasta_file);


/* call once at end of program */
void
chunk_strategy_free();


/* call if you are re-running a new chunk of input that is within a
  thread_queue_run() call. resets everything except n_all_loci_read,
  n_all_bytes_read. */
void
chunk_strategy_reset();


/* call to define a new workload. sets span.  sets n_loci_total to the
   size of the intersection between span and the query_regions.  sets
   n_loci_read to zero. */
void
chunk_strategy_set_span(struct contig_span span);


/* estimate the bytes wanted based on the strategy */
unsigned long
cs_max_bytes_wanted(); 


#endif /* _CHUNK_STRATEGY_H */
