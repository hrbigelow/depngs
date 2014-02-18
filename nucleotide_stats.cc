#include "nucleotide_stats.h"
#include "stats_tools.h"
#include "pileup_tools.h"

#include <cstring>
#include <cassert>
#include <numeric>
#include <cmath>


NucleotideStats::NucleotideStats()
{
    size_t D = Nucleotide::num_bqs;

    this->jpd_buffer = new double[D * 4];
    this->cpd_buffer = new double[D * 4];

    for (size_t b = 0; b != 4; ++b)
    {
        this->complete_jpd[b] = this->jpd_buffer + (b * D);
        this->founder_base_likelihood[b] = this->cpd_buffer + (b * D);
    }
    // this->index_mapping = new std::string[D];
}

NucleotideStats::~NucleotideStats()
{
    delete this->jpd_buffer;
    delete this->cpd_buffer;
    // delete[] this->index_mapping;
}

namespace Nucleotide
{
    int const base_to_index[] =
        {
            4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,
            4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,
            4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,
            4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,
            4,0,4,1,4,4,4,2,4,4,4,4,4,4,4,4,
            4,4,4,4,3,4,4,4,4,4,4,4,4,4,4,4,
            4,0,4,1,4,4,4,2,4,4,4,4,4,4,4,4,
            4,4,4,4,3,4,4,4,4,4,4,4,4,4,4,4,
            4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,
            4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,
            4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,
            4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,
            4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,
            4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,
            4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,
            4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4
        };
    char const* bases_upper = "ACGTN";
    char const* strands = "+-";
    size_t const PLUS_STRAND = 0;
    size_t const MINUS_STRAND = 1;

    size_t const highest_quality = 50;
    size_t const num_s = 2;
    size_t const num_qs = num_s * (highest_quality + 1);
    size_t const num_bqs = 4 * num_qs;
    
    size_t encode(char basecall, size_t quality, size_t strand_index)
    {
        size_t basecall_index = Nucleotide::base_to_index[static_cast<size_t>(basecall)];
        
        return 
            (basecall_index * num_qs) 
            + (quality * num_s) 
            + strand_index;
    }

    void decode(size_t code, char * basecall, size_t *quality, size_t *strand)
    {
        *basecall = Nucleotide::bases_upper[code / num_qs];
        *quality = (code % num_qs) / num_s; 
        *strand = (code % num_s) == 0 ? PLUS_STRAND : MINUS_STRAND;
    }


};


//initialize all distributions from the counts_map.  counts_map
//is not necessarily normalized.  any entries in this stats object
//without corresponding entries in counts_map are regarded as having
//zero counts

// fields initialized are:
// jpd_buffer, index_mapping, name_mapping, complete_jpd, founder_base_marginal,
// founder_base_likelihood
/*
void NucleotideStats::initialize(JPD_DATA const& counts_map)
{
    size_t D = Nucleotide::num_bqs;

    JPD_DATA::const_iterator cit;
    std::fill(this->jpd_buffer, this->jpd_buffer + (4 * D), 0.0);

    size_t datum_index;
    for (cit = counts_map.begin(), datum_index = 0; 
         cit != counts_map.end(); 
         ++cit, ++datum_index)
    {
        // this->index_mapping[datum_index] = (*cit).first;
        // this->name_mapping[(*cit).first] = datum_index;
        for (size_t b = 0; b != 4; ++b)
        {
            this->complete_jpd[b][datum_index] = (*cit).second.data[b];
        }
    }
    
    normalize(this->jpd_buffer, 4 * D, this->jpd_buffer);
    
    for (size_t b = 0; b != 4; ++b)
    {
        this->founder_base_marginal[b] =
            std::accumulate(this->complete_jpd[b],
                            this->complete_jpd[b] + D, 0.0);
    }
    
    for (size_t b = 0; b != 4; ++b)
    {
        for (size_t di = 0; di != D; ++di)
        {
            this->founder_base_likelihood[b][di] =
                this->complete_jpd[b][di]
                / this->founder_base_marginal[b];
        }
    }
    
}
*/


 /*
JPD_DATA parse_jpd_rdb_file(char const* rdb_file);


void NucleotideStats::initialize_from_file(char const* rdb_file)
{
    JPD_DATA raw_counts = parse_jpd_rdb_file(rdb_file);
    this->initialize(raw_counts);
}
 */



void NucleotideStats::initialize(char const* rdb_file)
{
    //intialize data_prior
    FILE * rdb_fh = fopen(rdb_file, "r");
    
    if (rdb_fh == NULL)
    {
        fprintf(stderr, "Couldn't open data jpd file %s\n", 
                rdb_file);
        exit(5);
    }

    size_t D = Nucleotide::num_bqs;
    std::fill(this->jpd_buffer, this->jpd_buffer + (4 * D), 0.0);

    // JPD_DATA counts_map;
    // char name[100];
    double counts[4];
    double counts_sum;


    char basecall;
    int quality;
    char strand;
    size_t index_code;

    while (! feof(rdb_fh))
    {
        fscanf(rdb_fh, "%c_%i_%c\t%lf\t%lf\t%lf\t%lf\n", &basecall, &quality, &strand, 
               counts, counts+1, counts+2, counts+3);

        for (size_t bi = 0; bi != 4; ++bi)
        {
            if (counts[bi] < 0)
            {
                fprintf(stderr, "NucleotideStats::parse_rdb_file: "
                        "found negative count for %c_%i_%c.\n", basecall, quality, strand);
                exit(11);
            }
        }
        counts_sum = counts[0] + counts[1] + counts[2] + counts[3];

        if (counts_sum == 0)
        {
            continue;
        }

        index_code = Nucleotide::encode(basecall, quality, 
                                             (strand == '+' ? Nucleotide::PLUS_STRAND
                                              : Nucleotide::MINUS_STRAND));

        for (size_t b = 0; b != 4; ++b)
        {
            this->complete_jpd[b][index_code] = counts[b];
        }
    }
    fclose(rdb_fh);


    normalize(this->jpd_buffer, 4 * D, this->jpd_buffer);
    
    for (size_t b = 0; b != 4; ++b)
    {
        this->founder_base_marginal[b] =
            std::accumulate(this->complete_jpd[b],
                            this->complete_jpd[b] + D, 0.0);
    }
    
    for (size_t b = 0; b != 4; ++b)
    {
        for (size_t di = 0; di != D; ++di)
        {
            this->founder_base_likelihood[b][di] =
                this->complete_jpd[b][di]
                / this->founder_base_marginal[b];
        }
    }
}


//creates a new NucleotideStats object with the same data set
//as the calling one, but with a marginal distribution of observed data
//reflected by 'locus', but the shape P(fb|obs) of the calling object
/*
JPD_DATA
NucleotideStats::make_per_locus_stats(PileupSummary const& locus)
{
    JPD_DATA counts_map;
    for (size_t di = 0; di != this->num_distinct_data; ++di)
    {
        double zero_counts[] = { 0, 0, 0, 0 };
        counts_map.insert(std::make_pair(this->index_mapping[di], nuc_frequency(zero_counts)));
    }

    double locus_slice[4];
    for (size_t raw_index = 0; raw_index != locus.num_distinct_data; ++raw_index)
    {
        size_t di = locus.stats_index[raw_index];
        double orig_marginal = 
            this->complete_jpd[0][di]
            + this->complete_jpd[1][di]
            + this->complete_jpd[2][di]
            + this->complete_jpd[3][di];

        double adjust_factor;
        if (orig_marginal == 0)
        {
            assert(locus.raw_counts[raw_index] == 0);
            adjust_factor = 1.0;
        }
        else
        {
            adjust_factor = locus.raw_counts[raw_index] / orig_marginal;
        }

        //double adjust_factor = 1.0;
        for (size_t bi = 0; bi != 4; ++bi)
        {
            locus_slice[bi] = this->complete_jpd[bi][di] * adjust_factor;
        }
        counts_map.erase(this->index_mapping[di]);
        counts_map.insert(std::make_pair(this->index_mapping[di], nuc_frequency(locus_slice)));
    }
    return counts_map;
}
*/


 // !!! check this.
 // the nuc_frequency array is indexed according to locus.stats_index
 /*
nuc_frequency *
NucleotideStats::per_locus_stats(packed_counts & locus)
{
    nuc_frequency * new_stats = new nuc_frequency[locus.num_distinct_data];

    for (size_t r = 0; r != pc.num_distinct_data; ++r)
    {
        size_t c = locus.stats_index[r];
        double marg =
            this->complete_jpd[0][c]
            + this->complete_jpd[1][c]
            + this->complete_jpd[2][c]
            + this->complete_jpd[3][c];
        double adjust_factor = (marg == 0) ? 1.0 : locus.raw_counts[r] / marg;

        for (size_t bi = 0; bi != 4; ++bi)
        {
            new_stats[c].data[bi] = this->complete_jpd[bi][di] * adjust_factor;
        }

    }
    return new_stats;
}
 */


// pack the subset of statistics for this model that are represented in
// the stats_index of packed_counts
void NucleotideStats::pack(packed_counts * c)
{
    double * buf = c->fbqs_cpd;
    size_t D = c->num_data;

    for (size_t i = 0; i != D; ++i)
    {
        for (size_t f = 0; f != 4; ++f)
        {
            (*buf) = this->founder_base_likelihood[f][i];
            ++buf;
        }
    }
}
