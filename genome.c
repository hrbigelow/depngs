#include "genome.h"

#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <limits.h>
#include <assert.h>

#include "htslib/faidx.h"
#include "ordering.h"
#include "cache.h"

struct refseq reference_seq;

/* initialize reference_seq */
void
genome_init(const char *fasta_file, unsigned do_load_seqs)
{
    faidx_t *fai = fai_load(fasta_file);
    if (fai == NULL) {
        fprintf(stderr, "Error: %s:%u: Couldn't open fasta index file %s.fai\n",
                __FILE__, __LINE__, fasta_file);
        exit(1);
    }

    reference_seq.n_contig = faidx_nseq(fai);
    reference_seq.names = kh_init(s);
    reference_seq.contig = 
        malloc(reference_seq.n_contig * sizeof(reference_seq.contig[0]));
    
    int len;
    unsigned t;
    khiter_t it;
    int ret;
    for (t = 0; t != reference_seq.n_contig; ++t) {
        strcpy(reference_seq.contig[t].name, faidx_iseq(fai, t));
        it = kh_put(s, reference_seq.names, reference_seq.contig[t].name, &ret);
        kh_val(reference_seq.names, it) = t;

        if (do_load_seqs) {
            reference_seq.contig[t].seq = 
                faidx_fetch_seq(fai, reference_seq.contig[t].name, 0, INT_MAX, &len);
            if (reference_seq.contig[t].seq == NULL) {
                fprintf(stderr, "Error: %s:%u: Didn't find sequence %s in fasta file %s\n",
                        __FILE__, __LINE__, reference_seq.contig[t].name, fasta_file);
                exit(1);
            }
        }
        else
            reference_seq.contig[t].seq = NULL;
    }
    fai_destroy(fai);
}


void
genome_free()
{
    unsigned t;
    for (t = 0; t != reference_seq.n_contig; ++t)
        if (reference_seq.contig[t].seq)
            free(reference_seq.contig[t].seq);
    free(reference_seq.contig);
    kh_destroy(s, reference_seq.names);
}


int
genome_contig_order(const char *contig)
{
    khiter_t it = kh_get(s, reference_seq.names, contig);
    if (it == kh_end(reference_seq.names))
        return -1;
    else
        return kh_val(reference_seq.names, it);
}



/* */
struct pair_ordering_range *
parse_query_ranges(const char *query_range_file,
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
        int tid = genome_contig_order(contig);
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
