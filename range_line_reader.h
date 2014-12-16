#ifndef _RANGE_LINE_READER_H
#define _RANGE_LINE_READER_H

#include "file_binary_search.h"

/* controls the behavior of the reader.  shared across multiple calls
   to the reader. */
struct range_line_reader_par {
    FILE *fh;
    struct file_bsearch_index *ix; /* the current index into the file */
    struct file_bsearch_range *q, *qend; /* collection of regions desired */
    get_line_ord_t get_line_ord;
    int new_query;
    size_t max_line_size;
};

void range_line_reader(void *par, char **in_buf, size_t *in_size, size_t *in_alloc);

#endif /* _RANGE_LINE_READER_H */
