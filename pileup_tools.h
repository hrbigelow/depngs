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

    PileupSummary(int indel_histo_size);
    ~PileupSummary();

    void load_line(char const* line);

    char reference[100];
    int position;
    char reference_base;
    size_t read_depth;
    // int indel_histo_size;
    int base_counts[num_base_symbols]; //ACGTNacgtn
    int base_qual_sums[num_base_symbols]; //qualities for corresponding counts
    int sum_of_counts;
    // int * indel_counts;
    char * bases;
    char * bases_upper;
    char * quality_codes;

    packed_counts counts;

    std::string const* index_mapping;

    std::map<std::string, int> * indel_seqs;

    size_t quality(size_t read_num) const;

    void parse(size_t min_quality_score);

    static void SetFtype(FastqType _fastq_type);
};

FastqType FastqFileType(char const* pileup_file,
                        char * chunk_buffer_in,
                        size_t chunk_size,
                        size_t num_threads);


#endif // _PILEUP_TOOLS_H
