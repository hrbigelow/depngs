#ifndef _LOCUS_H
#define _LOCUS_H

#include "file_binary_search.h"

struct pair_ordering init_locus(const char *line);

int less_locus_range(const void *pa, const void *pb);

struct pair_ordering_range *
parse_query_ranges(const char *query_range_file, size_t *num_queries);


#endif /* _LOCUS_H */
