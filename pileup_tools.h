#ifndef _PILEUP_TOOLS_H
#define _PILEUP_TOOLS_H

#include <map>
#include <cstring>
#include <string>

#include "tools.h"
#include "nucleotide_stats.h"
#include "cache.h"

/* A class for representing a single line of pileup information
 */
const int num_base_symbols = 10;

struct less_char_ptr
{
    bool operator()(const char *a, const char *b)
    {
        return strcmp(a, b) < 0;
    }
};

typedef std::map<const char*, unsigned, less_char_ptr> CHAR_MAP;

class PileupSummary {

public:

    static char code_to_redux[];
    static const char nucleotides[];
    static int quality_code_offset;
    static FastqType fastq_type;

    PileupSummary();
    ~PileupSummary();
    PileupSummary& operator=(const PileupSummary &);
    void load_line(char const* line);

    char reference[100], reference_base;
    int position;

    /* all reads that span this locus, including reads containing
       gaps, and low-quality bases */
    size_t read_depth; 

    /* the subset of 'read_depth' reads that contain a matching (CIGAR
       'M') base at this locus. */
    size_t read_depth_match;

    /* the subset of 'read_depth_match' reads that have above
       threshold quality score */
    size_t read_depth_high_qual; 

    int base_counts[num_base_symbols]; //ACGTNacgtn
    int sum_of_counts;
    struct managed_buf bases, bases_upper, bases_raw, quality_codes;

    packed_counts counts;
    CHAR_MAP insertions, deletions;

    size_t quality(size_t read_num) const;

    void parse(size_t min_quality_score);

    static void SetFtype(FastqType _fastq_type);
    static void set_offset(int offset);
};

FastqType FastqFileType(char const* pileup_file,
                        char * chunk_buffer_in,
                        size_t chunk_size,
                        size_t num_threads);


// return the probable fastq offset (33 or 64) for this pileup file.
// scan up to chunk_size bytes.
int fastq_offset(const char *pileup_file,
                 char *chunk_buffer_in,
                 size_t chunk_size);


#endif // _PILEUP_TOOLS_H
