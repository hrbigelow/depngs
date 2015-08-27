#include "file_utils.h"

extern "C" {
#include "common_tools.h"
}


// return the probable fastq offset (33 or 64) for this pileup file.
// scan up to chunk_size bytes.
// 
int fastq_offset(const char *pileup_file,
                 char *chunk_buffer_in,
                 size_t chunk_size)
{
    unsigned char min = 255, max = 0;
    FILE *pileup_input_fh = open_if_present(pileup_file, "r");
    if (! pileup_input_fh)
        return -1;

    size_t fread_nsec;

    FileUtils::read_until_newline(chunk_buffer_in, chunk_size, 1e6, 
                                  pileup_input_fh, &fread_nsec);
    
    char *last_fragment;
    std::vector<char *> lines = 
        FileUtils::find_complete_lines_nullify(chunk_buffer_in, &last_fragment);

    std::vector<char *>::iterator it;
    PileupSummary locus;
    for (it = lines.begin(); it != lines.end(); ++it)
    {
        locus.load_line(*it);
        for (size_t rd = 0; rd != locus.read_depth_match; ++rd)
        {
            max = locus.quality_codes.buf[rd] > max ? locus.quality_codes.buf[rd] : max;
            min = locus.quality_codes.buf[rd] < min ? locus.quality_codes.buf[rd] : min;
        }
    }
    fclose(pileup_input_fh);

    // reject probable Solexa format
    if (min >= 59 && min < 64)
    {
        fprintf(stderr, "Warning: quality code minimum is %i, which is likely the unsupported Solexa format\n",
                min);
        return -1;
    }
    // Sanger and Illumina18 format range [33, 126]
    else if (min < 64 && max <= 126)
        return 33;

    // Illumina13 and Illumina15  have range [64, 126]
    else if (min >= 64 && max <= 126)
        return 64;

    else
    {
        fprintf(stderr, "Warning: quality code range is from %i to %i, which doesn't fit any known encoding\n",
                min, max);
        return -1;
    }

}
