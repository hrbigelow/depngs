#ifndef _PILEUP_TOOLS_H
#define _PILEUP_TOOLS_H

#include <map>
#include <string.h>
#include "tools.h"

extern "C" {
#include "nucleotide_stats.h"
#include "cache.h"
}


#define NUM_BASE_SYMBOLS 5

/* A class for representing a single line of pileup information
 */
struct less_char_ptr
{
    bool operator()(const char *a, const char *b)
    {
        return strcmp(a, b) < 0;
    }
};

typedef std::map<const char*, unsigned, less_char_ptr> CHAR_MAP;

void pileup_init(unsigned _min_qual_score);

class PileupSummary {

public:

    static const char nucleotides[];
    static int quality_code_offset;
    static FastqType fastq_type;

    PileupSummary();
    ~PileupSummary();

    void load_line(char *line);

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

    int base_counts_high_qual[NUM_BASE_SYMBOLS]; //ACGTN
    struct managed_buf bases, bases_upper, bases_raw, quality_codes;

    packed_counts counts;
    CHAR_MAP insertions, deletions;

    size_t quality(size_t read_num) const;

    void parse();

    static void SetFtype(FastqType _fastq_type);
    static void set_offset(int offset);
};

FastqType FastqFileType(const char *pileup_file,
                        char * chunk_buffer_in,
                        size_t chunk_size,
                        size_t num_threads);


// return the probable fastq offset (33 or 64) for this pileup file.
// scan up to chunk_size bytes.
int fastq_offset(const char *pileup_file,
                 char *chunk_buffer_in,
                 size_t chunk_size);


#endif // _PILEUP_TOOLS_H
