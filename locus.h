#ifndef _LOCUS_H
#define _LOCUS_H

#include "ordering.h"

/* useful as a function for getting an ordinal from a pileup
   line. must call genome_init(fasta_file, 0) before using. */
struct pair_ordering
init_locus(const char *line);

int
less_locus_range(const void *pa, const void *pb);

#endif /* _LOCUS_H */
