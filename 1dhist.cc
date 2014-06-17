//  parse a pileup file, record a depth histogram of each chromosome

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <vector>
#include <cassert>
#include <math.h>

#include "samutil/file_utils.h"

// return the bin number that the <val> falls into
// bin 0 is [0, 1/B), center is 1/2B
// bin 1 is [1/B, 2/B), center is 3/2B
// etc
size_t bin(float val, size_t bins_per_unit)
{
    return (size_t) floor(val * bins_per_unit);
}

// return the real coordinate of the center of the bin of a given size
float bin_center(size_t bin, size_t bins_per_unit)
{
    return ((bin * 2l) + 1l) / (2 * (float)bins_per_unit);
}

int main(int argc, char **argv)
{
    FILE * fh = fopen(argv[1], "r");
    if (! fh)
    {
        fprintf(stderr, "error, couldn't open input file %s\n", argv[1]);
        exit(1);
    }
    size_t bufsize = static_cast<size_t>(atof(argv[2]));
    size_t bins_per_unit = static_cast<size_t>(atof(argv[3]));
    size_t nunits = static_cast<size_t>(atof(argv[4]));
    size_t nbins = bins_per_unit * nunits;

    char *tag = argv[5];

    size_t max_pileup_line_size = 1e7;
    // char *chunk = new char[bufsize + 1];
    char *chunk = (char *)malloc(sizeof(char) * (bufsize + max_pileup_line_size + 1));
    char contig[10] = "";
    char prev_contig[10] = "";
    float depth;

    size_t *hist = (size_t *)calloc(nbins, sizeof(size_t));
    //size_t *hist = new size_t[nbins];
    size_t fread_nsec, bytes_read;
    char *last_fragment;
    int ssval;
    
    size_t bytes_wanted = bufsize;
    size_t dbin;

    while (!feof(fh))
    {
        bytes_read = FileUtils::read_until_newline(chunk, bytes_wanted, max_pileup_line_size, fh, &fread_nsec);
        std::vector<char *> lines_vec = FileUtils::find_complete_lines_nullify(chunk, &last_fragment);
        size_t nlines = lines_vec.size();
        char **lines = lines_vec.data();
        char **line;
        for (line = lines; line != lines + nlines; ++line)
        {
            ssval = sscanf(*line, "chr%s\t%*u\t%*c\t%f", contig, &depth);
            assert(ssval == 2);
            dbin = bin(depth, bins_per_unit);

            if (strcmp(prev_contig, contig) && strlen(prev_contig) != 0)
            {
                // print histogram
                for (size_t bin = 0; bin != nbins; ++bin)
                {
                    printf("%s\t%s\t%5.3f\t%Zu\n", tag, prev_contig, bin_center(bin, bins_per_unit), hist[bin]);
                    hist[bin] = 0;
                }
                fflush(stdout);
            }
            if (dbin < nbins)
            {
                hist[dbin]++;
            }
            strcpy(prev_contig, contig);
        }
    }

    // print out last chromosome histogram
    for (size_t bin = 0; bin != nbins; ++bin)
    {
        printf("%s\t%s\t%5.3f\t%Zu\n", tag, prev_contig, bin_center(bin, bins_per_unit), hist[bin]);
        hist[bin] = 0;
    }

    // delete[] chunk;
    // delete[] hist;
    free(chunk);
    free(hist);
    fclose(fh);
}
