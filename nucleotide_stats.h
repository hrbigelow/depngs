#ifndef _NUCLEOTIDE_STATS_H
#define _NUCLEOTIDE_STATS_H

#define NUC_HIGHEST_QUALITY 94
#define NUC_NUM_S 2
#define NUC_NUM_QS (2 * (NUC_HIGHEST_QUALITY + 1))
#define NUC_NUM_BQS (4 * NUC_NUM_QS)
#define NUC_NUM_FBQS (4 * NUC_NUM_BQS)
#define NUC_PLUS_STRAND 0
#define NUC_MINUS_STRAND 1

#include <stdlib.h>


/* holds the overall statistics for one sample. */
struct nucleotide_stats {
    double jpd[4][NUC_NUM_BQS]; /* P(F,B,Q,S).  Ordered F,B,Q,S */
    double cpd[4][NUC_NUM_BQS]; /* P(B,Q,S|F). */
    double marg[4]; /* ?? Sum over BQS? */
};

struct cpd_count {
    /* holds P(b,q,s|f) for each f, and particular (b,q,s) */
    double cpd[4];
    unsigned long ct; /* the count of this (b,q,s) tuple */
};

/* this structure holds a summary of all raw base calls at a given
   locus. the values in stats_index can be decoded into (basecall,
   quality, strand) triplets using Nucleotide::decode */
struct packed_counts
{
    struct cpd_count stats[NUC_NUM_BQS];
    /* number of populated elements in stats and 'stats_index'. */
    size_t num_data; 

    /* stats_index[n] = code.  decode_nucleotide(code, &b, &q,
       &s). gives the (b,q,s) tuple. */
    size_t stats_index[NUC_NUM_BQS];
};

/* initialize default phred error probabilities in global variable
   'error_probability' */
void
init_phred();

int base_to_index(char base);
void nucleotide_stats_initialize(const char *rdb_file, struct nucleotide_stats *s);
extern size_t encode_nucleotide(char basecall, size_t quality, size_t strand);
extern void decode_nucleotide(size_t code, char *basecall, size_t *quality, size_t *strand);
void nucleotide_stats_pack(const struct nucleotide_stats *stats, struct packed_counts *c);

#endif // _NUCLEOTIDE_STATS_H
