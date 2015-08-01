#ifndef _LOCUS_H
#define _LOCUS_H

#include "ordering.h"


/* initialize fasta index resources.  call at start of program. */
void
locus_init(const char *fasta_file);

/* free fasta index resources.  call at end of program. */
void
locus_free();

/* useful as a function for getting an ordinal from a pileup
   line. must call genome_init(fasta_file, 0) before using. */
struct pair_ordering
parse_pileup_locus(const char *line);

#endif /* _LOCUS_H */
