#ifndef _NUCLEOTIDE_STATS_H
#define _NUCLEOTIDE_STATS_H

#include <map>
#include <string>

namespace Nucleotide 
{
    //transforms ACGT and acgt to 0123, everything else to 4
    extern int const base_to_index[];
    extern const char *bases_upper;
    extern const char *strands;
    extern const size_t highest_quality;
    extern const size_t num_s;
    extern const size_t num_qs;
    extern const size_t num_bqs;
    extern const size_t PLUS_STRAND;
    extern const size_t MINUS_STRAND;

    extern size_t encode(char basecall, size_t quality, size_t strand);
    extern void decode(size_t code, char *basecall, size_t *quality, size_t *strand);

};

// this structure holds a summary of all raw base calls at a given locus
// the values in stats_index can be decoded into (basecall, quality, strand) triplets
// using Nucleotide::decode
struct packed_counts
{
    unsigned long *raw_counts;
    size_t *stats_index;
    double *fbqs_cpd; // founder base likelihood in f,b,q,s order
    size_t num_data; // number of elements in 'raw_counts' and 'stats_index'.  fbqs_cpd has 4 * raw_counts elements.
};


class PileupSummary;


/* 
 */
class NucleotideStats {

 public:
    // P(founder_base, basecall, quality, strand).  Order is F,B,Q,S
    double *jpd_buffer;

    // P(basecall, quality, strand | founder_base).  Order is F,B,Q,S
    double *cpd_buffer;

    double founder_base_marginal[4];
    double *complete_jpd[4];
    double *founder_base_likelihood[4];

    NucleotideStats();
    ~NucleotideStats();
    void initialize(const char *rdb_file);

    void pack(packed_counts *c);
};


#endif // _NUCLEOTIDE_STATS_H
