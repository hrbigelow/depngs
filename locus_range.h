#ifndef _LOCUS_RANGE_H
#define _LOCUS_RANGE_H

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
struct contig_region {
    unsigned tid;
    unsigned beg, end;
};


struct contig_pos {
    unsigned tid;
    unsigned pos;
};


struct contig_span {
    struct contig_pos beg, end;
};


/* 0 if there is any overlap.  -1 if a is less, 1 if a is greater */
int
cmp_contig_region(const struct contig_region a, 
                  const struct contig_region b);


int
cmp_contig_pos(const struct contig_pos a,
               const struct contig_pos b);


struct contig_region *
parse_locus_ranges(const char *locus_range_file,
                   unsigned *n_queries,
                   unsigned long *n_total_loci);


#endif /* _LOCUS_RANGE_H */
