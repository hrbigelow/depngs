#include "locus.h"
#include "defs.h"
#include "fasta.h"

#include <assert.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>


/* initialize fasta index resources */
void
locus_init(const char *fasta_file)
{
    fasta_thread_init(fasta_file);
}


void
locus_free()
{
    fasta_thread_free();
}

/* initialize a locus from a character line */
struct pair_ordering
parse_pileup_locus(const char *line)
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
    int tid;
    if ((tid = fasta_get_tid(contig)) >= 0)
        o.hi = (size_t)tid;
    else {
        fprintf(stderr, "%s:%u: Error: Couldn't find contig %s in fasta index\n",
                __FILE__, __LINE__, contig);
        exit(1);
    }

    o.lo = (size_t)pos;
    return o;
}
