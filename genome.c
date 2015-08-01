#include "genome.h"

#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <limits.h>
#include <assert.h>

#include "htslib/faidx.h"
#include "cache.h"

struct refseq reference_seq;

/* initialize reference_seq */
void
genome_init(const char *fasta_file)
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
genome_load_contig(unsigned tid)
{
    if (tid >= reference_seq.n_contig) {
        fprintf(stderr, "Error: tid %u must be less than the number of contigs (%u)\n",
                tid, reference_seq.n_contig);
        return 1;
    }
    if (reference_seq.contig[tid].seq != NULL) {
        fprintf(stderr, "Error: reference sequence %s (tid %u) already loaded\n",
                reference_seq.contig[tid].name, tid);
        return 1;
    }
    reference_seq.contig[tid].seq = 
        faidx_fetch_seq(fai, reference_seq.contig[tid].name, 0, INT_MAX, &len);

    if (reference_seq.contig[tid].seq == NULL) {
        fprintf(stderr, "Error: %s:%u: Didn't find sequence %s in fasta file %s\n",
                __FILE__, __LINE__, reference_seq.contig[tid].name, fasta_file);
        return 1;
    }
    return 0;
}


int
genome_free_contig(unsigned tid)
{
    if (tid >= reference_seq.n_contig) {
        fprintf(stderr, "Error: tid %u must be less than the number of contigs (%u)\n",
                tid, reference_seq.n_contig);
        return 1;
    }
    if (reference_seq.contig[tid].seq == NULL) {
        fprintf(stderr, "Error: reference sequence %s (tid %u) not loaded (or already freed)\n",
                reference_seq.contig[tid].name, tid);
        return 1;
    }
    free(reference_seq.contig[tid].seq);
    reference_seq.contig[tid].seq = NULL;
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
