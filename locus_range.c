#include "locus_range.h"
#include "fasta.h"
#include "cache.h"

#include <stdlib.h>
#include <limits.h>
#include <stdio.h>
#include <assert.h>


#define CMP(a, b) ((a) == (b) ? 0 : ((a) < (b) ? -1 : 1))

#define MIN(a, b) ((a) < (b) ? (a) : (b))
#define MAX(a, b) ((a) < (b) ? (b) : (a))


int
cmp_contig_region(const struct contig_region a, 
                  const struct contig_region b)
{
    int cmp;
    if ((cmp = CMP(a.tid, b.tid)) != 0) return cmp;
    if (CMP(a.end, b.beg) < 1) return -1;
    if (CMP(b.end, a.beg) < 1) return 1;
    return 0;
}


int
cmp_contig_pos(const struct contig_pos a,
               const struct contig_pos b)
{
    int cmp;
    if ((cmp = CMP(a.tid, b.tid)) != 0) return cmp;
    return CMP(a.pos, b.pos);
}


static int
qsort_wrapper(const void *a, const void *b)
{
    const struct contig_region *ca = a, *cb = b;
    return cmp_contig_region(*ca, *cb);
}

struct contig_region *
parse_locus_ranges(const char *locus_range_file,
                   unsigned *n_queries,
                   unsigned long *n_total_loci)
{
    struct contig_region *queries;
    if (! locus_range_file) {
        *n_total_loci = 0;
        unsigned c, n = fasta_nseq();
        queries = malloc(n * sizeof(struct contig_region));
        *n_queries = n;
        int len;
        for (c = 0; c != n; ++c) {
            len = fasta_seq_ilen(c);
            queries[c] = (struct contig_region){ c, 0, len };
            *n_total_loci += len;
        }
        return queries;
    }
        
    FILE *locus_range_fh = fopen(locus_range_file, "r");
    if (! locus_range_fh) {
        fprintf(stderr, "Error: %s: Couldn't open %s\n",
                __func__, locus_range_file);
        exit(1);
    }
        
    size_t n_alloc = 10;
    queries = malloc(n_alloc * sizeof(struct contig_region));
        
    /* construct the set of non-overlapping query ranges */
    char contig[1000];
    unsigned beg_pos, end_pos, rval, nq = 0;

    while (!feof(locus_range_fh)) {
        if ((rval = fscanf(locus_range_fh, "%s\t%u\t%u\n",
                           contig, &beg_pos, &end_pos)) != 3) {
            fprintf(stderr, "Error: %s: bad format in file %s\n",
                    __func__, locus_range_file);
            exit(1);
        }
        int tid = fasta_get_tid(contig);
        if (tid == -1) {
            fprintf(stderr, "%s:%u: Error: Couldn't find contig %s in fasta index\n",
                    __FILE__, __LINE__, contig);
            exit(1);
        }
        queries[nq++] = (struct contig_region){ tid, beg_pos, end_pos };
        ALLOC_GROW(queries, nq + 1, n_alloc);
    }   
    fclose(locus_range_fh);

    /* non-overlapping intervals are well ordered.  overlaps compare
       equal and thus the ordering is undefined. */
    qsort(queries, nq, sizeof(queries[0]), qsort_wrapper);
    
    /* resolve overlaps. p = previous, q = current, w = write pointer.
       NOTE: K&R p 103: "The address of the first element past the end
       of an array can be used in pointer arithmetic".  Therefore the
       loop test 'q < queries + nq' will properly fail if nq == 0 and q
       == queries + 1) */
    struct contig_region *q, *p = NULL, *w;
    for (p = w = queries, q = p + 1; q < queries + nq; ++q) {
        if (cmp_contig_region(*p, *q) == 0) {
            p->beg = MIN(p->beg, q->beg);
            p->end = MAX(p->end, q->end);
        } else {
            *w++ = *p;
            p = q;
        }
    }
    *n_queries = w - queries;
    for (q = queries, *n_total_loci = 0; q != queries + *n_queries; ++q)
        *n_total_loci += q->end - q->beg;

    return queries;
}
