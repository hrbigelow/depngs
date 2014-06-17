/*
  Output the first four fields of a pileup file, also stripping the
  'chr' string from the contig
 */


#include <vector>
#include <string.h>

#include "samutil/file_utils.h"

int main(int argc, char *argv[])
{
    if (argc != 4)
    {
        printf("%s <bufsize> <input_file> <output_file>\n", argv[0]);
        exit(1);
    }
    size_t bufsize = static_cast<size_t>(atof(argv[1]));
    FILE *in_fh = fopen(argv[2], "r");
    FILE *out_fh = fopen(argv[3], "w");

    if (! in_fh)
    {
        fprintf(stderr, "Error: couldn't open input file %s\n", argv[2]);
        exit(1);
    }
    if (! out_fh)
    {
        fprintf(stderr, "Error: couldn't open output file %s\n", argv[3]);
        exit(1);
    }

    size_t bytes_wanted = bufsize, bytes_read, max_pileup_line_size = 1e7, fread_nsec;
    char *chunk = (char *)malloc(sizeof(char) * (bufsize + max_pileup_line_size + 1));
    char *last_fragment;
    // size_t position, depth;
    // char chrom[100];
    // char refbase;
    char *end;

    while (! feof(in_fh))
    {
        bytes_read = FileUtils::read_until_newline(chunk, bytes_wanted, max_pileup_line_size, in_fh, &fread_nsec);
        std::vector<char *> lines_vec = FileUtils::find_complete_lines_nullify(chunk, &last_fragment);
        size_t nlines = lines_vec.size();
        char **lines = lines_vec.data();
        char **line;
        
        for (line = lines; line != lines + nlines; ++line)
        {
            *line += 3; // get rid of 'chr'
            end = strchr(strchr(strchr(strchr(*line, '\t') + 1, '\t') + 1, '\t') + 1, '\t');
            *end++ = '\n';
            fwrite(*line, 1, end - *line, out_fh);
        }
    }
    
    delete[] chunk;
    fclose(in_fh);
    fclose(out_fh);
}
