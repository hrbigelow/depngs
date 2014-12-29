#include "pileup_bsearch.h"
#include "dict.h"

#include <assert.h>
#include <stdlib.h>

#define MISSING_CONTIG(CTG)                         \
    do {                                            \
    fprintf(stderr, "Error at %s: %u. "             \
            "Contig %s not found in dictionary.\n", \
            __FILE__, __LINE__, (contig));          \
    exit(1);                                        \
} while (0)
      
      
/* initialize a locus from a character line */
struct pair_ordering init_locus(const char *line)
{
    char contig[200];
    unsigned pos;
    struct pair_ordering o;
    int nparsed = sscanf(line, "%s\t%u\t", contig, &pos);
    assert(nparsed == 2);
    long ix;
    if ((ix = dict_search(contig)) >= 0)
        o.hi = (size_t)ix;
    else
        MISSING_CONTIG(contig);
    o.lo = (size_t)pos;
    return o;
}
