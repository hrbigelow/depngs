#include "chunk_strategy.h"

#include <stdlib.h>

struct chunk_strategy cs_stats;

#define MAX(a,b) ((a) < (b) ? (b) : (a))

/* in order to simplify the interface, we choose a reasonable default
   value for this. it may be off by a factor of 10 or so, but that
   won't matter for large input. */
#define DEFAULT_BYTES_PER_LOCUS 100

void
chunk_strategy_init(unsigned n_files, unsigned n_threads)
{
    cs_stats.pos = (struct contig_pos){ 0, 0 };
    cs_stats.n_files = n_files;
    cs_stats.n_threads = n_threads;
    cs_stats.n_bytes_read = calloc(n_files, sizeof(cs_stats.n_bytes_read[0]));
    cs_stats.n_loci_read = 0;
    cs_stats.default_bytes_per_locus = DEFAULT_BYTES_PER_LOCUS;
}


void
chunk_strategy_free()
{
    free(cs_stats.n_bytes_read);
}


void
chunk_strategy_reset(unsigned long n_loci)
{
    cs_stats.pos = (struct contig_pos){ 0, 0 };
    cs_stats.n_loci_total = n_loci;
    cs_stats.n_loci_read = 0;
    unsigned s;
    for (s = 0; s != cs_stats.n_files; ++s)
        cs_stats.n_bytes_read[s] = 0;
}


void
cs_stats_reset_pos()
{
    cs_stats.pos = (struct contig_pos){ 0, 0 };
}


/* void cs_set_defaults(unsigned long max_bytes_small_chunk, */
/*                      unsigned long small_chunk_size, */
/*                      unsigned long default_bytes_per_locus) */
/* { */
/*     cs_stats.max_bytes_small_chunk = max_bytes_small_chunk; */
/*     cs_stats.small_chunk_size = small_chunk_size; */
/*     cs_stats.default_bytes_per_locus = default_bytes_per_locus; */
/* } */

/* return the estimated maximum number of bytes that we should attempt
   to parse per thread. */
unsigned long
cs_max_bytes_wanted()
{
    /* find maximum of n_bytes_per_locus across samples */
    unsigned long max_bytes_per_locus = 0;
    if (cs_stats.n_loci_read == 0)
        max_bytes_per_locus = cs_stats.default_bytes_per_locus;
    else {
        unsigned long max_bytes_read = 0;
        unsigned f;
        for (f = 0; f != cs_stats.n_files; ++f)
            max_bytes_read = MAX(max_bytes_read, cs_stats.n_bytes_read[f]);
        max_bytes_per_locus = max_bytes_read / cs_stats.n_loci_read;
    }

    /* the most bytes that are left in any one sample */
    unsigned long most_bytes_left = 
        (cs_stats.n_loci_total - cs_stats.n_loci_read)
        * max_bytes_per_locus;

    return most_bytes_left / cs_stats.n_threads;
}


unsigned cs_bytes_per_locus(unsigned s)
{
    return cs_stats.n_bytes_read[s] / cs_stats.n_loci_read;
}
