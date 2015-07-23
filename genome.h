#ifndef _GENOME_H
#define _GENOME_H

#include "htslib/faidx.h"
#include "khash.h"

#define REFNAME_MAXLEN 300

KHASH_MAP_INIT_STR(s, unsigned);

struct refseq {
    struct {
        char name[REFNAME_MAXLEN + 1];
        char *seq;
    } *contig;
    unsigned n_contig;
    khash_t(s) *names;
};

extern struct refseq reference_seq;

/* initialize reference_seq */
void
genome_init(const char *fasta_file, unsigned do_load_seqs);

void
genome_free();

/* return the ordering of the contig, or -1 if not found */
int
genome_contig_order(const char *contig);

struct pair_ordering_range *
parse_query_ranges(const char *query_range_file,
                   unsigned *num_queries,
                   unsigned long *n_total_loci);


#endif /* _GENOME_H */
