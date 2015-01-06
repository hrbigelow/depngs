#ifndef _NUCLEOTIDE_STATS_H
#define _NUCLEOTIDE_STATS_H

#include <map>
#include <string>

#define NUC_HIGHEST_QUALITY 94
#define NUC_NUM_S 2
#define NUC_NUM_QS (2 * (NUC_HIGHEST_QUALITY + 1))
#define NUC_NUM_BQS (4 * NUC_NUM_QS)
#define NUC_NUM_FBQS (4 * NUC_NUM_BQS)

namespace Nucleotide 
{
    //transforms ACGT and acgt to 0123, everything else to 4
    extern int const base_to_index[];
    extern const char *bases_upper;
    extern const char *strands;
    /* extern const size_t highest_quality; */
    /* extern const size_t num_s; */
    /* extern const size_t num_qs; */
    /* extern const size_t num_bqs; */
    extern const size_t PLUS_STRAND;
    extern const size_t MINUS_STRAND;

    extern size_t encode(char basecall, size_t quality, size_t strand);
    extern void decode(size_t code, char *basecall, size_t *quality, size_t *strand);

};

/* this structure holds a summary of all raw base calls at a given
   locus the values in stats_index can be decoded into (basecall,
   quality, strand) triplets using Nucleotide::decode */
struct packed_counts
{
    /* raw_counts[n] = number of occurrences of a particular (b,q,s)
       tuple at this locus */
    unsigned long *raw_counts;

    /* stats_index[n] = code.  use Nucleotide::decode(code, &b, &q,
       &s). use this together with raw_counts[n] to find the number of
       occurrences of (b,q,s) tuples. */
    size_t *stats_index; 

    /* fbqs_cpd[nf] gives P(b,q,s|f).  nf / 4 corresponds to n.  nf %
       4 corresponds to founder base (A,C,G,T). */
    double *fbqs_cpd;

    /* number of elements in 'raw_counts' and 'stats_index'.  fbqs_cpd
       has 4 * raw_counts elements. */
    size_t num_data; 
};


class PileupSummary;


/* 
 */
class NucleotideStats {

 public:
    // P(founder_base, basecall, quality, strand).  Order is F,B,Q,S
    double jpd_buffer[NUC_NUM_FBQS];

    // P(basecall, quality, strand | founder_base).  Order is F,B,Q,S
    double cpd_buffer[NUC_NUM_FBQS];

    double founder_base_marginal[4];
    double *complete_jpd[4];
    double *founder_base_likelihood[4];

    NucleotideStats();
    ~NucleotideStats();
    void initialize(const char *rdb_file);

    void pack(packed_counts *c);
};


#endif // _NUCLEOTIDE_STATS_H
