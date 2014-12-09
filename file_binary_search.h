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

#include <unistd.h>
#include <stdio.h>

/* a generic space in which to store line ordinal information */
struct file_bsearch_ord {
    size_t hi, lo;
};


/* allows mapping file offsets to loci */
struct file_bsearch_index {
    struct { 
        struct file_bsearch_ord beg, end; 
    } span;
    off_t start_offset, end_offset;
    struct file_bsearch_index *left, *right, *parent;
    char *span_contents; /* if non-null, will contain contents of
                            file */
};

/* initialize static function pointers */
void file_bsearch_init(struct file_bsearch_ord (*_get_line_ord)(const char *line),
                       size_t _mem_scan_threshold);


/* defines the ordering logic for the file_bsearch_ord structure. */
int less_file_bsearch_ord(const void *pa, const void *pb);


/* generate a root index that spans the whole file */
struct file_bsearch_index *find_root_index(FILE *fh);


/* cuts the index range in two, finding the nearest linebreak after
   the offset midpoint, and then choosing the portion that contains
   cur. uses file seeking until a predefined size, then reads the
   entire range into span_contents.  stops when the nearest linebreak
   is at the end of the range. */
struct file_bsearch_index *
find_loose_index(struct file_bsearch_index *ix, struct file_bsearch_ord cur, FILE *fh);


/* release resources.  call this after you are done searching, or
   before you want to call file_bsearch_init for a different kind of
   data.  */
void file_bsearch_free();


/* return offset position of start of line where cur is (or would be
   inserted) */
off_t off_lower_bound(const struct file_bsearch_index *ix, struct file_bsearch_ord cur);

/* return offset position of start of line just after cur (or just
   after where cur would be inserted) */
off_t off_upper_bound(const struct file_bsearch_index *ix, struct file_bsearch_ord cur);


/* free the index tree */
void file_bsearch_index_free(struct file_bsearch_index *root);


#endif /* _FILE_BINARY_SEARCH_H */
