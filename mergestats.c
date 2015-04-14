#include <stdio.h>
#include <getopt.h>
#include <stdlib.h>
#include <sys/stat.h>
#include <string.h>

#include "metropolis_sampling.h"
#include "dirichlet_diff_cache.h"

int mergestats_usage()
{
    fprintf(stderr,
            "\nUsage: dep mergestats [options] diststats1.rdb diststats2.rdb ...\n"
            "Options:\n\n"
            "-o STRING     OUTFILE, merged diststats file\n"
            "\n"
            "Merges several input diststats files, checking that the headers match\n"
            "\n"
            );
    return 1;
}

#define MIN(a,b) ((a) < (b) ? (a) : (b))

void cat(char *buf, 
         size_t buf_size, 
         size_t bytes_to_print, 
         FILE *in_fh, 
         FILE *out_fh)
{
    while (! feof(in_fh))
    {
        size_t bytes_chunk = MIN(bytes_to_print, buf_size);
        fread(buf, 1, bytes_chunk, in_fh);
        fwrite(buf, 1, bytes_chunk, out_fh);
        bytes_to_print -= bytes_chunk;
        if (bytes_to_print == 0) break;
    }
    fflush(out_fh);
}


int main_mergestats(int argc, char **argv)
{
    char *out_file = NULL;
    char c;
    while ((c = getopt(argc, argv, "o:")) >= 0)
    {
        switch(c)
        {
        case 'o': out_file = optarg; break;
        default: return mergestats_usage(); break;
        }
    }
    size_t n_infiles = argc - optind;
    if (n_infiles < 2) return mergestats_usage();

    FILE *out_fh = fopen(out_file, "w");
    if (! out_fh)
    {
        fprintf(stderr, "Error: Couldn't open output file %s for writing\n", out_file);
        exit(1);
    }

    struct posterior_settings pset, psetn;
    memset(&pset, 0, sizeof(pset));
    memset(&psetn, 0, sizeof(psetn));

    unsigned max1, max2, max1n, max2n;
    size_t copy_buf_size = 1e8;
    char *copy_buf = malloc(copy_buf_size);

    struct stat s;
    size_t file_size, hdr_size;
    FILE *in_fh;

    unsigned i;
    for (i = 0; i != n_infiles; ++i)
    {
        in_fh = fopen(argv[i + optind], "r");
        if (! in_fh)
        {
            fprintf(stderr, "Error: Couldn't open input file %s for reading\n", argv[i]);
            exit(1);
        }
        file_size = (fstat(fileno(in_fh), &s), s.st_size);
    
        if (i == 0)
        {
            parse_diststats_header(in_fh, &pset, &max1, &max2);
            rewind(in_fh);
            cat(copy_buf, copy_buf_size, file_size, in_fh, out_fh);
        }
        else
        {
            parse_diststats_header(in_fh, &psetn, &max1n, &max2n);
            if (memcmp(&pset, &psetn, sizeof(pset)) != 0
                || max1 != max1n
                || max2 != max2n)
            {
                fprintf(stderr, "Error: input file %u has a different header than first input file\n", i + 1);
                exit(1);
            }
            hdr_size = (size_t)ftell(in_fh);
            cat(copy_buf, copy_buf_size, file_size - hdr_size, in_fh, out_fh);
        }
        fclose(in_fh);
    }
    fclose(out_fh);
    free(copy_buf);

    return 0;
}
