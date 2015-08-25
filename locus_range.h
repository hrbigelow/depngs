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


#define CONTIG_REGION_BEG(r) (struct contig_pos){ (r).tid, (r).beg }
#define CONTIG_REGION_END(r) (struct contig_pos){ (r).tid, (r).end }

#define MIN_CONTIG_POS(a, b) (cmp_contig_pos((a), (b)) < 0 ? (a) : (b))
#define MAX_CONTIG_POS(a, b) (cmp_contig_pos((a), (b)) < 0 ? (b) : (a))

/* 0 if there is any overlap.  -1 if a is less, 1 if a is greater */
int
cmp_contig_region(const struct contig_region a, 
                  const struct contig_region b);


int
cmp_contig_pos(const struct contig_pos a,
               const struct contig_pos b);


/* parse locus range file in the format of 'contig <tab> start <tab>
   end'.  assume start and end are given in 1-based coordinates. store
   as zero-based coordinates. */
struct contig_region *
parse_locus_ranges(const char *locus_range_file,
                   const char *fasta_file,
                   unsigned *n_queries,
                   unsigned long *n_total_loci);


/* find the subrange of the sorted range [qbeg, qend) that intersects
   subset, storing it in *qlo, *qhi. return the total number of loci
   in the intersection. */
unsigned long
find_intersecting_span(const struct contig_region *qbeg,
                       const struct contig_region *qend,
                       struct contig_span subset,
                       const struct contig_region **qlo,
                       const struct contig_region **qhi);


/* return the largest contig position end such that the intersection
   of [beg, end) and [qbeg, qend) is <= n_max_loci. sets *n_found_loci
   to the number of loci that are actually in the intersection of
   [qbeg, qend) and [beg, end) */
struct contig_pos
find_span_of_size(const struct contig_region *qbeg,
                  const struct contig_region *qend,
                  struct contig_pos beg,
                  unsigned long n_max_loci,
                  unsigned long *n_found_loci);


/* returns the contig_region that is the intersection between r and s,
   or a zero-length region if there is no intersection. */
struct contig_region
region_span_intersect(struct contig_region r, struct contig_span s);


#endif /* _LOCUS_RANGE_H */
