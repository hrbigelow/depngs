#ifndef _PILEUP_BSEARCH_H
#define _PILEUP_BSEARCH_H

/* adapter functions for doing binary search on pileup files, via the
   generic file_binary_search mechanism. */

#include "file_binary_search.h"

struct file_bsearch_ord init_locus(const char *line);

struct locus_range {
    struct file_bsearch_ord beg, end;
};

int less_locus_range(const void *pa, const void *pb);          

#endif /* _PILEUP_BSEARCH_H */