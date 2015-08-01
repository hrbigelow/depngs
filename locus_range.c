#include "locus_range.h"
#include "fasta.h"
#include "cache.h"

#include <stdlib.h>
#include <limits.h>
#include <stdio.h>
#include <assert.h>

struct pair_ordering_range *
parse_locus_ranges(const char *locus_range_file,
                   unsigned *n_queries,
                   unsigned long *n_total_loci)
{
    struct pair_ordering_range *queries;
    if (! locus_range_file) {
        queries = malloc(sizeof(struct pair_ordering_range));
        queries[0] = 
            (struct pair_ordering_range){ 
            { 0, 0 }, { SIZE_MAX - 2, 0 } 
        };
        *n_queries = 1;
        *n_total_loci = ULONG_MAX;
        return queries;
    }
        
    FILE *locus_range_fh = fopen(locus_range_file, "r");
    if (! locus_range_fh) {
        fprintf(stderr, "Error: %s: Couldn't open %s\n",
                __func__, locus_range_file);
        exit(1);
    }
        
    size_t n_alloc = 10;
    queries = malloc(n_alloc * sizeof(struct pair_ordering_range));
        
    /* construct the set of non-overlapping query ranges */
    char contig[1000];
    unsigned beg_pos, end_pos, rval;
    *n_queries = 0;
    *n_total_loci = 0;

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
        
        queries[*n_queries] = (struct pair_ordering_range) {
            { tid, beg_pos }, { tid, end_pos }
        };
            
        ++*n_queries;
        ALLOC_GROW_TYPED(queries, *n_queries + 1, n_alloc);
    }   
    fclose(locus_range_fh);

    qsort(queries, *n_queries, sizeof(queries[0]), cmp_pair_ordering_range);
    
    struct pair_ordering_range *q, *p = NULL;
    for (q = queries; q != queries + *n_queries - 1; ++q) {
        if (p && cmp_pair_ordering(&p->end, &q->beg) > 0) {
            /* must be on same contig if they are overlapping and
               sorted */
            assert(p->end.hi == q->beg.hi);
            q->beg.lo = p->end.lo;
            if (q->end.lo < q->beg.lo)
                q->end.lo = q->beg.lo;
        }

        p = q;
    }
    for (q = queries; q != queries + *n_queries; ++q)
        *n_total_loci += q->end.lo - q->beg.lo;

    return queries;
}
