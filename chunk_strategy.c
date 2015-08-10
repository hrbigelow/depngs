#include "chunk_strategy.h"

#include <stdlib.h>

struct chunk_strategy cs_stats;

#define MAX(a,b) ((a) < (b) ? (b) : (a))

/* call this if a range file is given */
void cs_init_by_range(unsigned n_loci_total, unsigned n_files)
{
    cs_stats.pos = (struct contig_pos){ 0, 0 };
    cs_stats.do_range_estimation = 1;
    cs_stats.n_bytes_total = NULL;
    cs_stats.n_bytes_read = calloc(n_files, sizeof(cs_stats.n_bytes_read[0]));
    cs_stats.n_loci_total = n_loci_total;
    cs_stats.n_loci_read = 0;
}


/* call this if no range file is given */
void cs_init_whole_file(unsigned n_files)
{
    cs_stats.pos = (struct contig_pos){ 0, 0 };
    cs_stats.do_range_estimation = 0;
    cs_stats.n_bytes_total = malloc(n_files * sizeof(cs_stats.n_bytes_total[0]));
    cs_stats.n_bytes_read = calloc(n_files, sizeof(cs_stats.n_bytes_read[0]));
    cs_stats.n_loci_total = 0;
    cs_stats.n_loci_read = 0;
}


void cs_set_total_bytes(unsigned i, unsigned long bytes)
{
    cs_stats.n_bytes_total[i] = bytes;
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
   see .h file for notes.
 */
unsigned cs_get_bytes_wanted(unsigned n_files)
{
    unsigned f;
    unsigned long most_bytes_left = 0;
    if (cs_stats.do_range_estimation) {
        for (f = 0; f != n_files; ++f) {
            unsigned n_bytes_per_locus = 
                cs_stats.n_loci_read == 0
                ? cs_stats.default_bytes_per_locus
                : cs_stats.n_bytes_read[f] / cs_stats.n_loci_read;
            
            most_bytes_left = MAX((cs_stats.n_loci_total - cs_stats.n_loci_read)
                                  * n_bytes_per_locus, 
                                  most_bytes_left);
        }
    }
    else
        for (f = 0; f != n_files; ++f)
            most_bytes_left = MAX(most_bytes_left,
                                  cs_stats.n_bytes_total[f]
                                  - cs_stats.n_bytes_read[f]);

    unsigned bytes_wanted = 
        most_bytes_left > cs_stats.max_bytes_small_chunk
        ? most_bytes_left
        : cs_stats.small_chunk_size;

    return bytes_wanted;
}


unsigned cs_bytes_per_locus(unsigned s)
{
    return cs_stats.n_bytes_read[s] / cs_stats.n_loci_read;
}

void cs_free()
{
    if (cs_stats.n_bytes_total) free(cs_stats.n_bytes_total);
    free(cs_stats.n_bytes_read);
}
