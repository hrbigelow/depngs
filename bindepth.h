#ifndef _BINDEPTH_H
#define _BINDEPTH_H

/*

bindepth format is:

file: contig-dictionary contig-data [, contig-data ...]

  contig-dictionary: num-contigs contig [, contig]
    num-contigs: size_t representing the number of contigs
    contig: contig-name contig-size
    junk: uninitialized bytes filling up to dict-bytes

      contig-name: 20 char, null-terminated character string
      contig-size: size_t representing contig size

   contig-data: depth [, depth]
      depth: float representing the depth.

 */

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
