#ifndef _PILEUP_TOOLS_H
#define _PILEUP_TOOLS_H

#include <map>
#include <cstring>
#include <string>

#include "tools.h"
#include "nucleotide_stats.h"

/* A class for representing a single line of pileup information
 */
const int num_base_symbols = 10;


class PileupSummary {

public:

    static char code_to_redux[];
    static char const nucleotides[];
    static int quality_code_offset;
    static FastqType fastq_type;

 PileupSummary(int indel_histo_size) : indel_histo_size(indel_histo_size),
        bases(NULL), bases_upper(NULL), quality_codes(NULL)
    {
        counts.raw_counts = NULL;
        counts.stats_index = NULL;
        counts.fbqs_cpd = NULL;
        memset(base_counts, 0, sizeof(base_counts[0]) * num_base_symbols);
        memset(base_qual_sums, 0, sizeof(base_counts[0]) * num_base_symbols);
        indel_counts = new int[indel_histo_size * 2 + 1];
        indel_seqs = new std::map<std::string, int>[indel_histo_size * 2 + 1];
        indel_counts += indel_histo_size;
        indel_seqs += indel_histo_size;
        for (int i = -indel_histo_size; i <= indel_histo_size; ++i)
        {
            indel_counts[i] = 0;
            indel_seqs[i] = std::map<std::string, int>();
        }
        sum_of_counts = 0;
    }


    ~PileupSummary()
    {
        delete &indel_counts[-indel_histo_size];
        delete [] &indel_seqs[-indel_histo_size];
        if (bases != NULL)
        {
            delete bases;
            bases = NULL;
        }
        if (bases_upper != NULL)
        {
            delete bases_upper;
            bases_upper = NULL;
        }
        if (quality_codes != NULL)
        {
            delete quality_codes;
            quality_codes = NULL;
        }
        if (this->counts.raw_counts != NULL)
        {
            delete this->counts.raw_counts;
            this->counts.raw_counts = NULL;
        }
        if (this->counts.stats_index != NULL)
        {
            delete this->counts.stats_index;
            this->counts.stats_index = NULL;
        }
    }

    void load_line(char const* line);

    char reference[100];
    int position;
    char reference_base;
    size_t read_depth;
    int indel_histo_size;
    int base_counts[num_base_symbols]; //ACGTNacgtn
    int base_qual_sums[num_base_symbols]; //qualities for corresponding counts
    int sum_of_counts;
    int * indel_counts;
    char * bases;
    char * bases_upper;
    char * quality_codes;

    packed_counts counts;

    std::string const* index_mapping;

    std::map<std::string, int> * indel_seqs;

    size_t quality(size_t read_num) const;

    void parse(size_t min_quality_score);

    FastqType FastqFileType(char const* pileup_file,
                            char * chunk_buffer_in,
                            size_t chunk_size);

    static void SetFtype(FastqType _fastq_type);
};

#endif // _PILEUP_TOOLS_H
