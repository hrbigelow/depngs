#include "chunk_strategy.h"

#include <stdlib.h>

struct chunk_strategy cs_stats;

#define MAX(a,b) ((a) < (b) ? (b) : (a))

void
chunk_strategy_init(unsigned n_files)
{
    cs_stats.pos = (struct contig_pos){ 0, 0 };
    cs_stats.n_files = n_files;
    cs_stats.n_bytes_read = calloc(n_files, sizeof(cs_stats.n_bytes_read[0]));
    cs_stats.n_loci_read = 0;
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


void cs_set_defaults(unsigned long max_bytes_small_chunk,
                     unsigned long small_chunk_size,
                     unsigned long default_bytes_per_locus)
{
    cs_stats.max_bytes_small_chunk = max_bytes_small_chunk;
    cs_stats.small_chunk_size = small_chunk_size;
    cs_stats.default_bytes_per_locus = default_bytes_per_locus;
}

/* estimate the bytes wanted based on the strategy.
   see .h file for notes. */
unsigned long
cs_get_bytes_wanted(unsigned n_files)
{
    unsigned long most_bytes_left = 0;
    unsigned f;
    for (f = 0; f != n_files; ++f) {
        unsigned n_bytes_per_locus = 
            cs_stats.n_loci_read == 0
            ? cs_stats.default_bytes_per_locus
            : cs_stats.n_bytes_read[f] / cs_stats.n_loci_read;
        
        most_bytes_left = MAX((cs_stats.n_loci_total - cs_stats.n_loci_read)
                              * n_bytes_per_locus, 
                              most_bytes_left);
    }
    unsigned long bytes_wanted = 
        most_bytes_left > cs_stats.max_bytes_small_chunk
        ? most_bytes_left
        : cs_stats.small_chunk_size;

    return bytes_wanted;
}


unsigned cs_bytes_per_locus(unsigned s)
{
    return cs_stats.n_bytes_read[s] / cs_stats.n_loci_read;
}
