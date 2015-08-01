#ifndef _LOCUS_RANGE_H
#define _LOCUS_RANGE_H

#include "ordering.h"

/* 
   before using this function, fasta_init(fasta_file) must be called,
   using the fasta_file relevant to the set of locus ranges.

   parse a locus range file.  translates contig names to tids.
   merges overlapping ranges, and orders them.

   locus_range_file:  has lines of <contig>\t<start>\t<end>

   fasta_file: the matching faidx index file <fasta_file>.fai must be
   present.  It is used to initialize the fasta index which allows
   mapping contig names to tids.

   returns a pair_ordering_range which must be freed by the caller.
*/
struct pair_ordering_range *
parse_locus_ranges(const char *locus_range_file,
                   unsigned *n_queries,
                   unsigned long *n_total_loci);


#endif /* _LOCUS_RANGE_H */
