#ifndef _BINDEPTH_H
#define _BINDEPTH_H

#include <cstddef>
#include <cstdio>

struct contig_dict_t
{
    char name[10];
    char shortname[10];
    size_t size;
};

contig_dict_t *parse_contig_dict(FILE *contig_dict_fh, size_t *ncontigs);
void write_contig_dict(contig_dict_t *contigs, size_t ncontigs, FILE *out_fh);
contig_dict_t *read_contig_dict(FILE *in_fh, size_t *ncontigs);

#endif // _BINDEPTH_H
