/* Functions for binary search on a sorted, line-oriented text file.
   The client provides a function for extracting a file_bsearch_ord
   structure from any given line in the file.  The client then
   provides a query file_bsearch_ord structure representing a desired
   line position, and the library then uses this extraction function
   to do binary search on the file (using seeking) to return the
   physical offset (off_t) where the given query would reside.  (The
   actual line doesn't have to exist in the file).

   struct file_bsearch_ord represents the index of the given line, and
   is ordered as:

   a.hi < b.hi || (a.hi == b.hi && a.lo < b.lo)
*/

#ifndef _FILE_BINARY_SEARCH_H
#define _FILE_BINARY_SEARCH_H

#ifndef _GNU_SOURCE
#define _GNU_SOURCE /* needed for memrchr */
#endif

#define _FILE_OFFSET_BITS 64

#include <unistd.h>
#include <stdio.h>

#include "ordering.h"

/* allows mapping file offsets to loci */
struct file_bsearch_node {
    struct pair_ordering span_beg, span_end; 
    off_t start_offset, end_offset;
    struct file_bsearch_node *left, *right, *parent;
    char *span_contents; /* if non-null, will contain contents of
                            file */
};

/* represents the index tree 'as a whole', together with the
   associated structure(s) needed to use it. */
struct file_bsearch_index {
    FILE *fh;
    char *mem_scan_buf;
    char *line_buf;
    size_t line_len;
    struct file_bsearch_node *root, *cur_node;
    size_t n_nodes;
};

typedef struct pair_ordering (*get_line_ord_t)(const char *line);

/* initialize static function pointers */
void file_bsearch_init(get_line_ord_t _get_line_ord,
                       size_t _mem_scan_threshold);


/* generate an index that spans the whole file */
struct file_bsearch_index file_bsearch_make_index(const char *file);


/* release resources.  call this after you are done searching, or
   before you want to call file_bsearch_init for a different kind of
   data.  */
void file_bsearch_free();


/* return an estimated upper size limit to contain the logical range
   [beg, end) */
size_t range_to_size(struct file_bsearch_index *ix,
                     struct pair_ordering beg,
                     struct pair_ordering end);


/* return an estimated highegst upper bound logical position that uses
   <= size bytes */
struct pair_ordering size_to_range(struct file_bsearch_index *ix,
                                   struct pair_ordering beg,
                                   size_t size);


/* read a logical range into a buffer */
size_t read_range(struct file_bsearch_index *ix,
                  struct pair_ordering beg,
                  struct pair_ordering end,
                  char *buf);

/* return the largest offset such that all lines spanning query start
   at or after this offset. */
off_t off_lower_bound(struct file_bsearch_index *ix,
                      struct pair_ordering query);

/* free a portion of the index tree that covers [beg, end) range. */
void file_bsearch_index_free(struct file_bsearch_index ix);

/* free all nodes that are contained in [beg, end)*/
size_t file_bsearch_range_free(struct file_bsearch_index *ix,
                               struct pair_ordering beg,
                               struct pair_ordering end);


#endif /* _FILE_BINARY_SEARCH_H */
