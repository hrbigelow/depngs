#include "locus.h"
#include "cache.h"
#include "defs.h"
#include "htslib/faidx.h"

#include <assert.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>


static faidx_t *global_fai;

void
locus_global_init(const char *fai_file)
{
    global_fai = fai_load(fai_file);
}

void
locus_global_free()
{
    fai_destroy(global_fai);
}


/* initialize a locus from a character line */
struct pair_ordering
init_locus(const char *line)
{
    char contig[200];
    unsigned pos;
    struct pair_ordering o;

    char *p;
    unsigned ref_end = (p = (char *)memchr(line, '\t', MAX_PILEUP_LINE)) - line;
    assert(ref_end < sizeof(contig));
    strncpy(contig, line, ref_end);
    contig[ref_end] = '\0';

    pos = strtol(++p, &p, 10);
    assert(*p == '\t');

    /* int nparsed = sscanf(line, "%s\t%u\t", contig, &pos); */
    /* assert(nparsed == 2); */
    long ix;
    if ((ix = faidx_tid(contig)) >= 0)
        o.hi = (size_t)ix;
    else {
        fprintf(stderr, "%s:%u: Error: Couldn't find contig %s in fasta index\n",
                __FILE__, __LINE__, contig);
        exit(1);
    }

    o.lo = (size_t)pos;
    return o;
}



/* */
struct pair_ordering_range *
parse_query_ranges(const char *query_range_file,
                   const faidx_t *fai,
                   unsigned *num_queries,
                   unsigned long *n_total_loci)
{
    FILE *query_range_fh = fopen(query_range_file, "r");
    if (! query_range_fh) {
        fprintf(stderr, "Error: %s: Couldn't open %s\n",
                __func__, query_range_file);
        exit(1);
    }
        
    size_t num_alloc = 10;
    struct pair_ordering_range *queries = 
        (struct pair_ordering_range *)
        malloc(num_alloc * sizeof(struct pair_ordering_range));
        
    /* construct the set of non-overlapping query ranges */
    char contig[1000];
    unsigned beg_pos, end_pos, rval;
    *num_queries = 0;
    *n_total_loci = 0;

    while (!feof(query_range_fh)) {
        if ((rval = fscanf(query_range_fh, "%s\t%u\t%u\n",
                           contig, &beg_pos, &end_pos)) != 3) {
            fprintf(stderr, "Error: %s: bad format in file %s\n",
                    __func__, query_range_file);
            exit(1);
        }
        int tid = faidx_tid(fai, contig);
        if (tid == -1) {
            fprintf(stderr, "%s:%u: Error: Couldn't find contig %s in fasta index\n",
                    __FILE__, __LINE__, contig);
            exit(1);
        }
        
        queries[*num_queries] = (struct pair_ordering_range) {
            { tid, beg_pos }, { tid, end_pos }
        };
            
        ++*num_queries;
        ALLOC_GROW_TYPED(queries, *num_queries + 1, num_alloc);
    }   
    fclose(query_range_fh);

    qsort(queries, *num_queries, sizeof(queries[0]), cmp_pair_ordering_range);
    
    struct pair_ordering_range *q, *p = NULL;
    for (q = queries; q != queries + *num_queries - 1; ++q) {
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
    for (q = queries; q != queries + *num_queries; ++q)
        *n_total_loci += q->end.lo - q->beg.lo;

    return queries;
}
