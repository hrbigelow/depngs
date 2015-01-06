#include "nucleotide_stats.h"
#include "stats_tools.h"
#include "pileup_tools.h"

#include <cstring>
#include <cassert>
#include <numeric>
#include <cmath>


NucleotideStats::NucleotideStats()
{
    size_t D = NUC_NUM_BQS;
    for (size_t b = 0; b != 4; ++b)
    {
        this->complete_jpd[b] = this->jpd_buffer + (b * D);
        this->founder_base_likelihood[b] = this->cpd_buffer + (b * D);
    }
    // this->index_mapping = new std::string[D];
}

NucleotideStats::~NucleotideStats()
{
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
    const char *bases_upper = "ACGTN";
    const char *strands = "+-";
    const size_t PLUS_STRAND = 0;
    const size_t MINUS_STRAND = 1;

    // const size_t highest_quality = 94;
    // const size_t num_s = 2;
    // const size_t num_qs = num_s * (highest_quality + 1);
    // const size_t num_bqs = 4 * num_qs;
    
    size_t encode(char basecall, size_t quality, size_t strand_index)
    {
        size_t basecall_index = Nucleotide::base_to_index[static_cast<size_t>(basecall)];
        assert(quality <= NUC_HIGHEST_QUALITY);

        return 
            (basecall_index * NUC_NUM_QS) 
            + (quality * NUC_NUM_S) 
            + strand_index;
    }

    void decode(size_t code, char * basecall, size_t *quality, size_t *strand)
    {
        *basecall = Nucleotide::bases_upper[code / NUC_NUM_QS];
        *quality = (code % NUC_NUM_QS) / NUC_NUM_S; 
        *strand = (code % NUC_NUM_S) == 0 ? PLUS_STRAND : MINUS_STRAND;
    }


};




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
    std::fill(this->jpd_buffer, this->jpd_buffer + NUC_NUM_FBQS, 0.0);

    double counts[4], counts_sum;

    char basecall, strand;
    int quality;
    size_t index_code;

    while (! feof(rdb_fh))
    {
        fscanf(rdb_fh, "%c_%i_%c\t%lf\t%lf\t%lf\t%lf\n", &basecall, &quality, &strand, 
               counts, counts+1, counts+2, counts+3);

        for (size_t bi = 0; bi != 4; ++bi)
            if (counts[bi] < 0)
            {
                fprintf(stderr, "NucleotideStats::parse_rdb_file: "
                        "found negative count for %c_%i_%c.\n", basecall, quality, strand);
                exit(11);
            }

        counts_sum = counts[0] + counts[1] + counts[2] + counts[3];

        if (counts_sum == 0)
            continue;

        index_code = Nucleotide::encode(basecall, quality, 
                                             (strand == '+' ? Nucleotide::PLUS_STRAND
                                              : Nucleotide::MINUS_STRAND));

        for (size_t b = 0; b != 4; ++b)
            this->complete_jpd[b][index_code] = counts[b];
    }
    fclose(rdb_fh);

    normalize(this->jpd_buffer, NUC_NUM_FBQS, this->jpd_buffer);
    
    for (size_t b = 0; b != 4; ++b)
    {
        this->founder_base_marginal[b] =
            std::accumulate(this->complete_jpd[b],
                            this->complete_jpd[b] + NUC_NUM_BQS, 0.0);
    }
    
    for (size_t b = 0; b != 4; ++b)
        for (size_t di = 0; di != NUC_NUM_BQS; ++di)
            this->founder_base_likelihood[b][di] =
                this->complete_jpd[b][di]
                / this->founder_base_marginal[b];
}



// initialize fbqs_cpd field of c.  packs the subset of statistics for
// this model that are represented in the c->stats_index
void NucleotideStats::pack(packed_counts *c)
{
    double *buf = c->fbqs_cpd;
    size_t D = c->num_data;
    size_t code;
    for (size_t i = 0; i != D; ++i)
    {
        code = c->stats_index[i];
        for (size_t f = 0; f != 4; ++f)
        {
            *buf = this->founder_base_likelihood[f][code];
            assert(*buf != 0.0);
            ++buf;
        }
    }
}
