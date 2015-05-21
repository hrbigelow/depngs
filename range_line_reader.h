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

    /* total number of loci to process. if we don't have individual
       ranges, set this to ULONG_MAX.*/
    unsigned long n_loci_left; 
    struct pair_ordering start_pos, end_pos;
    struct off_pair *offset_pairs;
    unsigned int n_offsets;
    get_line_ord_t get_line_ord;
    int new_query;
};

/* read data starting at beg, using ranges and file handles specified
   in par.  store retrieved data in bufs */
void rl_reader(void *par, struct managed_buf *bufs);

/* set beg to the next read beginning position, without actually
   reading data or modifying bufs */
void rl_scanner(void *par, unsigned max_bytes);


struct reader_state {
    struct pair_ordering pos;
    unsigned long n_loci_left;
};

/* getter and setter for thread_queue's global_reader_state
   parameter. */
void rl_set_global_state(void *par, void *state);

void rl_get_global_state(void *par, void *state);


#endif /* _RANGE_LINE_READER_H */
