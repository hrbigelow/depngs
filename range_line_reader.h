#ifndef _RANGE_LINE_READER_H
#define _RANGE_LINE_READER_H

#include "file_binary_search.h"
#include "cache.h"

/* controls the behavior of the reader.  shared across multiple calls
   to the reader. */
struct range_line_reader_par {
    struct file_bsearch_index *ix; /* the current indexes into the files */
    size_t n_ix;

    /* pointers into the underlying shared range array */
    struct pair_ordering_range *qbeg, *qend;

    struct off_pair *offset_pairs;
    unsigned int n_offsets;
    get_line_ord_t get_line_ord;
    int new_query;
};

/* read data starting at beg, using ranges and file handles specified
   in par.  store retrieved data in bufs */
void rl_reader(void *par, struct managed_buf *bufs);

/* update global statistics and position markers to prepare for next
   read. this function only called by a single thread at a time */
void rl_scanner(void *par, unsigned max_bytes);

/* defaults for the chunking strategy */
/* If we have less than 1GB of input to go, switch to using small
   chunks of SMALL_CHUNK size */
#define MAX_BYTES_SMALL_CHUNK 1e9

/* This size is big enough for a pileup line and small enough to be a
   quickly-processed chunk, increasing the chance that multiple
   threads can keep working until the end of the input. */
#define SMALL_CHUNK 5e6

/* If in do_range_estimation mode, and we don't yet have an estimate
   for n_bytes_per_locus, use this as a default value. */
#define DEFAULT_BYTES_PER_LOCUS 100


#endif /* _RANGE_LINE_READER_H */
