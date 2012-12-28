#ifndef _PILEUP_TOOLS_H
#define _PILEUP_TOOLS_H

#include <map>
#include <cstring>
#include <string>

#include "tools.h"

/* A class for representing a single line of pileup information
 */
const int num_base_symbols = 10;


class PileupSummary {

public:

    static char code_to_redux[];
    static char const nucleotides[];
    static int quality_code_offset;
    static FastqType fastq_type;

 PileupSummary(int indel_histo_size) : _indel_histo_size(indel_histo_size),
        _bases(NULL), _bases_upper(NULL), _quality_codes(NULL)
    {
        memset(base_counts, 0, sizeof(base_counts[0]) * num_base_symbols);
        memset(base_qual_sums, 0, sizeof(base_counts[0]) * num_base_symbols);
        indel_counts = new int[indel_histo_size * 2 + 1];
        indel_seqs = new std::map<std::string, int>[indel_histo_size * 2 + 1];
        indel_counts += _indel_histo_size;
        indel_seqs += _indel_histo_size;
        for (int i = -indel_histo_size; i <= indel_histo_size; ++i)
        {
            indel_counts[i] = 0;
            indel_seqs[i] = std::map<std::string, int>();
        }
        sum_of_counts = 0;
    }
    ~PileupSummary()
    {
        delete &indel_counts[-_indel_histo_size];
        delete [] &indel_seqs[-_indel_histo_size];
        if (_bases != NULL)
        {
            delete _bases;
            _bases = NULL;
        }
        if (_bases_upper != NULL)
        {
            delete _bases_upper;
            _bases_upper = NULL;
        }
        if (_quality_codes != NULL)
        {
            delete _quality_codes;
            _quality_codes = NULL;
        }
    }

    bool load_line(FILE * file);

    char _reference[100];
    int _position;
    char _reference_base;
    size_t _read_depth;
    int _indel_histo_size;
    int base_counts[num_base_symbols]; //ACGTNacgtn
    int base_qual_sums[num_base_symbols]; //qualities for corresponding counts
    int sum_of_counts;
    int * indel_counts;
    char * _bases;
    char * _bases_upper;
    char * _quality_codes;

    std::map<std::string, int> * indel_seqs;

    size_t quality(size_t read_num) const;

    FastqType FastqFileType(char const* pileup_file);
    static void SetFtype(FastqType _fastq_type);
};

#endif // _PILEUP_TOOLS_H
