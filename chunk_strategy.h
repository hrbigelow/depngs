#ifndef _CHUNK_STRATEGY_H
#define _CHUNK_STRATEGY_H

#include "ordering.h"

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
    struct pair_ordering pos;

    unsigned do_range_estimation;

    unsigned long *n_bytes_total;
    unsigned long *n_bytes_read;
    unsigned long n_loci_total;
    unsigned long n_loci_read;

    unsigned long max_bytes_small_chunk; /* # bytes left to switch to
                                              small chunk strategy */
    unsigned long small_chunk_size;
    unsigned long default_bytes_per_locus;
};

/* call this if a range file is given */
void cs_init_by_range(unsigned n_loci_total, unsigned n_files);

/* call this if no range file is given */
void cs_init_whole_file(unsigned n_files);

/* call this for each file whose total bytes is known */
void cs_set_total_bytes(unsigned i, unsigned long bytes);

/* must call this */
void cs_set_defaults(unsigned long max_bytes_small_chunk,
                     unsigned long small_chunk_size,
                     unsigned long default_bytes_per_locus);


/* estimate the bytes wanted based on the strategy */
unsigned cs_get_bytes_wanted(unsigned n_files); 

void cs_free();

#endif /* _CHUNK_STRATEGY_H */
