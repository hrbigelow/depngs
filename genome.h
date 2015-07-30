#ifndef _GENOME_H
#define _GENOME_H

#include "htslib/faidx.h"
#include "khash.h"

#define REFNAME_MAXLEN 300

#define REFERENCE_SAMPLE UINT_MAX

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


/* parse query_range_file.  if NULL, return a single maximal range,
   setting n_queries to 1 and n_total_loci to a maximal value */
struct pair_ordering_range *
parse_query_ranges(const char *query_range_file,
                   unsigned *n_queries,
                   unsigned long *n_total_loci);


#endif /* _GENOME_H */
