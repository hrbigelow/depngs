#include "bindepth.h"

#include <cstdlib>
#include <string.h>

contig_dict_t *parse_contig_dict(FILE *contig_dict_fh, size_t *ncontigs)
{
    size_t space = 5;
    contig_dict_t *contigs = (contig_dict_t *)malloc(sizeof(contig_dict_t) * space);
    contig_dict_t *cdict = contigs;
    *ncontigs = 0;

    while (! feof(contig_dict_fh))
    {
        if ((size_t)(cdict - contigs) == space)
        {
            space *= 2;
            contigs = (contig_dict_t *)realloc(contigs, sizeof(contig_dict_t) * space);
            cdict = contigs + *ncontigs;
        }
        fscanf(contig_dict_fh, "%s\t%zu\n", cdict->name, &cdict->size);
        strcpy(cdict->shortname, cdict->name + (strcasestr(cdict->name, "chr") ? 3 : 0));
        ++cdict;
        ++(*ncontigs);
    }
    return contigs;
}


void write_contig_dict(contig_dict_t *contigs, size_t ncontigs, FILE *out_fh)
{
    fwrite(&ncontigs, sizeof(size_t), 1, out_fh);
    fwrite(contigs, sizeof(contigs[0]), ncontigs, out_fh);
}


// read a binary-encoded 
contig_dict_t *read_contig_dict(FILE *in_fh, size_t *ncontigs)
{
    fread(ncontigs, sizeof(size_t), 1, in_fh);
    contig_dict_t *contigs = (contig_dict_t *)malloc(sizeof(contig_dict_t) * *ncontigs);
    fread(contigs, sizeof(contig_dict_t), *ncontigs, in_fh);

    contig_dict_t *c = contigs, *cend = c + *ncontigs;
    while (c != cend)
    {
        strcpy(c->shortname, c->name + (strcasestr(c->name, "chr") ? 3 : 0));
        ++c;
    }
    return contigs;
}
