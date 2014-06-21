#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <string.h>

#include "bindepth.h"


struct locus_t
{
    char contig[10];
    size_t position, nbases;
};

int main(int argc, char **argv)
{
    if (argc != 4)
    {
        fprintf(stderr, "Usage: fasta_grep in.bfasta loci.rdb snippets.rdb\n");
        return 1;
    }

    const char *bfasta_file = argv[1], *loci_file = argv[2], *out_file = argv[3];

    FILE *bfasta_fh, *loci_fh, *out_fh;
    if (! (bfasta_fh = fopen(bfasta_file, "r")))
    {
        fprintf(stderr, "Error: couldn't open input file %s\n", bfasta_file);
        return 1;
    }
    if (! (loci_fh = fopen(loci_file, "r")))
    {
        fprintf(stderr, "Error: couldn't open input file %s\n", loci_file);
        return 1;
    }
    if (! (out_fh = fopen(out_file, "w")))
    {
        fprintf(stderr, "Error: couldn't open output file %s\n", out_file);
        return 1;
    }

    size_t ncontigs;
    contig_dict_t *contigs = read_contig_dict(bfasta_fh, &ncontigs), *contigp;

    size_t nloci = 10; // initial value
    locus_t *loci = (locus_t *)malloc(sizeof(locus_t) * nloci), *lp = loci, *lpe = lp + nloci;

    while (1)
    {
        int nparsed = fscanf(loci_fh, "%s\t%zu\t%zu\n", lp->contig, &lp->position, &lp->nbases);
        assert(nparsed == 3);
        ++lp;
        if (lp == lpe)
        {
            nloci *= 2;
            size_t lpos = lp - loci;
            loci = (locus_t *)realloc(loci, sizeof(locus_t) * nloci);
            lp = loci + lpos, lpe = loci + nloci;
        }
        int ch = fgetc(loci_fh);
        if (ch == EOF) { lpe = lp; lp = loci; break; }
        else ungetc(ch, loci_fh);
    }
    fclose(loci_fh);

    long datpos = ftell(bfasta_fh);
    char seq_buf[10000];
    while (lp != lpe)
    {
        // calculate the offset
        contigp = contigs;
        fseek(bfasta_fh, datpos, SEEK_SET); // reset to beginning of data
        while (strcmp(contigp->name, lp->contig)) { fseek(bfasta_fh, contigp++->size, SEEK_CUR); }
        fseek(bfasta_fh, lp->position, SEEK_CUR); // seek forward to desired position within contig
        fread(seq_buf, 1, lp->nbases, bfasta_fh);
        seq_buf[lp->nbases] = '\0';
        fprintf(out_fh, "%s\t%Zu\t%s\n", lp->contig, lp->position, seq_buf);
        ++lp;
    }
    
    fclose(bfasta_fh);
    fclose(out_fh);

    free((void *)contigs);
    free((void *)loci);

    return 0;
}
