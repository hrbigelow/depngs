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
    double jpd_buffer[NUC_NUM_FBQS]; /* P(F,B,Q,S).  Ordered F,B,Q,S */
    double cpd_buffer[NUC_NUM_FBQS]; /* P(B,Q,S|F). */

    double founder_base_marginal[4];
    double *complete_jpd[4];
    double *founder_base_likelihood[4];
};

struct cpd_count {
    double cpd[4]; /* = P(b,q,s|f) for each f in (A,C,G,T), and a
                      particular (b,q,s) */
    unsigned long ct; /* the count of this (b,q,s) tuple */
};

/* this structure holds a summary of all raw base calls at a given
   locus the values in stats_index can be decoded into (basecall,
   quality, strand) triplets using Nucleotide::decode */
struct packed_counts
{
    struct cpd_count stats[NUC_NUM_BQS];
    /* number of populated elements in stats and 'stats_index'. */
    size_t num_data; 

    /* stats_index[n] = code.  Nucleotide::decode(code, &b, &q,
       &s). gives the (b,q,s) tuple. */
    size_t stats_index[NUC_NUM_BQS];
};

int base_to_index(char base);
struct nucleotide_stats make_nucleotide_stats();
void nucleotide_stats_initialize(const char *rdb_file, struct nucleotide_stats *s);
extern size_t encode_nucleotide(char basecall, size_t quality, size_t strand);
extern void decode_nucleotide(size_t code, char *basecall, size_t *quality, size_t *strand);
void nucleotide_stats_pack(const struct nucleotide_stats *stats, struct packed_counts *c);

#endif // _NUCLEOTIDE_STATS_H
