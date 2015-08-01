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
genome_init(const char *fasta_file);

void
genome_free();

/* return the ordering of the contig, or -1 if not found */
int
genome_contig_index(const char *contig);

/* load a contig */
void
genome_load_contig(unsigned tid);


/* free a contig from memory */
void
genome_free_contig(unsigned tid);


#endif /* _GENOME_H */
