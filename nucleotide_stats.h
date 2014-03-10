#ifndef _NUCLEOTIDE_STATS_H
#define _NUCLEOTIDE_STATS_H

#include <map>
#include <string>

namespace Nucleotide 
{
    //transforms ACGT and acgt to 0123, everything else to 4
    extern int const base_to_index[];
    extern char const* bases_upper;
    extern char const* strands;
    extern size_t const highest_quality;
    extern size_t const num_s;
    extern size_t const num_qs;
    extern size_t const num_bqs;
    extern size_t const PLUS_STRAND;
    extern size_t const MINUS_STRAND;

    extern size_t encode(char basecall, size_t quality, size_t strand);
    extern void decode(size_t code, char * basecall, size_t *quality, size_t *strand);

};

// this structure holds a summary of all raw base calls at a given locus
// the values in stats_index can be decoded into (basecall, quality, strand) triplets
// using Nucleotide::decode
struct packed_counts
{
    unsigned long * raw_counts;
    size_t * stats_index;
    double * fbqs_cpd; // founder base likelihood in f,b,q,s order
    size_t num_data;
};


class PileupSummary;


/* 
 */
class NucleotideStats {

 public:
    // P(founder_base, basecall, quality, strand).  Order is F,B,Q,S
    double * jpd_buffer;

    // P(basecall, quality, strand | founder_base).  Order is F,B,Q,S
    double * cpd_buffer;

    double founder_base_marginal[4];
    double * complete_jpd[4];
    double * founder_base_likelihood[4];

    NucleotideStats();
    ~NucleotideStats();
    void initialize(char const* rdb_file);

    void pack(packed_counts * c);
};


#endif // _NUCLEOTIDE_STATS_H
