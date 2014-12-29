#ifndef _RANGE_LINE_READER_H
#define _RANGE_LINE_READER_H

#include "file_binary_search.h"
#include "cache.h"

/* controls the behavior of the reader.  shared across multiple calls
   to the reader. */
struct range_line_reader_par {
    struct file_bsearch_index *ix; /* the current indexes into the files */
    size_t n_ix;
    struct pair_ordering_range *q, *qend;
    get_line_ord_t get_line_ord;
    int new_query;
};



void range_line_reader(void *par, struct managed_buf *bufs);

#endif /* _RANGE_LINE_READER_H */
