#ifndef _RANGE_LINE_READER_H
#define _RANGE_LINE_READER_H

#include "file_binary_search.h"
#include "cache.h"

/* controls the behavior of the reader.  shared across multiple calls
   to the reader. */
struct range_line_reader_par {
    struct file_bsearch_index *ix; /* the current indexes into the files */
    size_t n_ix;
    struct pair_ordering_range *qcur, *qend;
    struct pair_ordering start_pos, end_pos;
    get_line_ord_t get_line_ord;
    int new_query;
};

/* read data starting at beg, using ranges and file handles specified
   in par.  store retrieved data in bufs */
void rl_reader(void *par, struct managed_buf *bufs);

/* set beg to the next read beginning position, without actually
   reading data or modifying bufs */
void rl_scanner(void *par, unsigned max_bytes);

/* set pos to the end_pos recorded in the par.  pos is a global value,
   so this function will be called within a mutex. */
void rl_set_start(void *par, void *pos);

void rl_get_start(void *par, void *pos);


#endif /* _RANGE_LINE_READER_H */
