#include "locus.h"
#include "defs.h"
#include "genome.h"

#include <assert.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>


/* initialize a locus from a character line */
struct pair_ordering
init_locus(const char *line)
{
    char contig[200];
    unsigned pos;
    struct pair_ordering o;

    char *p;
    unsigned ref_end = (p = (char *)memchr(line, '\t', MAX_PILEUP_LINE)) - line;
    assert(ref_end < sizeof(contig));
    strncpy(contig, line, ref_end);
    contig[ref_end] = '\0';

    pos = strtol(++p, &p, 10);
    assert(*p == '\t');

    /* int nparsed = sscanf(line, "%s\t%u\t", contig, &pos); */
    /* assert(nparsed == 2); */
    int ix;
    if ((ix = genome_contig_order(contig)) >= 0)
        o.hi = (size_t)ix;
    else {
        fprintf(stderr, "%s:%u: Error: Couldn't find contig %s in fasta index\n",
                __FILE__, __LINE__, contig);
        exit(1);
    }

    o.lo = (size_t)pos;
    return o;
}
