#ifndef _LOCUS_H
#define _LOCUS_H

#include "ordering.h"
#include "htslib/faidx.h"

/* call once at start of program */
void
locus_global_init(const char *fai_file);

/* call once at end of program to use */
void
locus_global_free();


/* useful as a function for getting an ordinal from a pileup line. */
struct pair_ordering
init_locus(const char *line);

int
less_locus_range(const void *pa, const void *pb);

struct pair_ordering_range *
parse_query_ranges(const char *query_range_file, 
                   const faidx_t *fai,
                   unsigned *num_queries,
                   unsigned long *n_total_loci);


#endif /* _LOCUS_H */
