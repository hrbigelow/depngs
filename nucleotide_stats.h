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

/*
struct nuc_frequency
{
    double data[4];
    nuc_frequency(double const _data[4])
    {
        std::copy(_data, _data + 4, this->data);
    }

    nuc_frequency(nuc_frequency const& a)
    {
        std::copy(a.data, a.data + 4, this->data);
    }

    double total() const
    {
        return data[0] + data[1] + data[2] + data[3];
    }
    void multiply(double factor)
    {
        data[0] *= factor;
        data[1] *= factor;
        data[2] *= factor;
        data[3] *= factor;
    }
};
*/    


// typedef std::map<std::string, nuc_frequency> JPD_DATA;

// this structure holds a summary of all raw base calls at a given locus
// the values in stats_index can be decoded into (basecall, quality, strand) triplets
// using Nucleotide::decode
struct packed_counts
{
    double * raw_counts;
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
    size_t num_distinct_data;
    /* std::map<std::string, size_t> name_mapping; */
    /* std::string * index_mapping; */

    NucleotideStats();
    ~NucleotideStats();
    /* void initialize(JPD_DATA const& counts_map); */
    void initialize(char const* rdb_file);
    // JPD_DATA make_per_locus_stats(PileupSummary const& locus);

    void pack(packed_counts * c);
};


#endif // _NUCLEOTIDE_STATS_H
