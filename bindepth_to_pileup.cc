// parse a bindepth file and output the truncated pileup format

#include <cstdlib>
#include "bindepth.h"

#define MIN(a,b) ((a) < (b) ? (a) : (b))

int main(int argc, char *argv[])
{
    if (argc == 1)
    {
        fprintf(stderr, "Usage: %s <bufsize> <infile.bindepth> <out.pileup>\n", argv[0]);
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

    size_t ncontigs;
    contig_dict_t *contigs = read_contig_dict(in_fh, &ncontigs);
    contig_dict_t *ctg;

    size_t pos, npos = bufsize / sizeof(float);
    float *buf = new float[npos];
    float *depth;
    for (ctg = contigs; ctg != contigs + ncontigs; ++ctg)
    {
        pos = 0;
        while (pos != ctg->size)
        {
            size_t n = MIN(npos, ctg->size - pos);
            fread(buf, sizeof(float), n, in_fh);
            depth = buf;
            for (; n != 0; ++pos, --n)
            {
                if (*depth != 0)
                {
                    fprintf(out_fh, "%s\t%zu\t%6.4f\n", ctg->name, pos + 1, *depth);
                }
                ++depth;
            }
        }
    }
    fclose(out_fh);
    fclose(in_fh);
    free(contigs);
    delete[] buf;
    return 0;
}
