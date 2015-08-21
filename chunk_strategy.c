#include "chunk_strategy.h"

#include <stdlib.h>

struct chunk_strategy cs_stats;

#define MAX(a,b) ((a) < (b) ? (b) : (a))

/* in order to simplify the interface, we choose a reasonable default
   value for this. it may be off by a factor of 10 or so, but that
   won't matter for large input. */
#define DEFAULT_BYTES_PER_LOCUS 100

void
chunk_strategy_init(unsigned n_files, unsigned n_threads,
                    unsigned long n_min_absolute_bytes,
                    const char *locus_range_file,
                    const char *fasta_file)
{
    cs_stats.n_files = n_files;
    cs_stats.n_threads = n_threads;
    cs_stats.n_min_absolute_bytes = n_min_absolute_bytes;
    cs_stats.n_all_bytes_read = calloc(n_files, sizeof(cs_stats.n_all_bytes_read[0]));
    cs_stats.n_all_loci_read = 0;
    cs_stats.query_regions = parse_locus_ranges(locus_range_file,
                                                fasta_file,
                                                &cs_stats.n_query_regions,
                                                &cs_stats.n_loci_total);
    /* initializes span and n_loci_total */
    chunk_strategy_reset();
}


void
chunk_strategy_free()
{
    free(cs_stats.n_all_bytes_read);
}


void
chunk_strategy_reset()
{
    struct contig_span span = 
        { CONTIG_REGION_BEG(cs_stats.query_regions[0]),
          CONTIG_REGION_END(cs_stats.query_regions[cs_stats.n_query_regions - 1]) };
    chunk_strategy_set_span(span);
}


void
chunk_strategy_set_span(struct contig_span span)
{
    cs_stats.span = span;
    const struct contig_region *qlo, *qhi;
    cs_stats.n_loci_total =
        find_intersecting_span(cs_stats.query_regions,
                               cs_stats.query_regions + cs_stats.n_query_regions,
                               span, &qlo, &qhi);
    cs_stats.n_loci_read = 0;
}


/* return the estimated maximum number of bytes that we should attempt
   to parse per thread. */
unsigned long
cs_max_bytes_wanted()
{
    /* find maximum of n_bytes_per_locus across samples */
    unsigned long max_bytes_per_locus = 0;
    if (cs_stats.n_all_loci_read == 0)
        max_bytes_per_locus = DEFAULT_BYTES_PER_LOCUS;
    else {
        unsigned long max_bytes_read = 0;
        unsigned f;
        for (f = 0; f != cs_stats.n_files; ++f)
            max_bytes_read = MAX(max_bytes_read, cs_stats.n_all_bytes_read[f]);
        max_bytes_per_locus = max_bytes_read / cs_stats.n_all_loci_read;
    }

    /* the most bytes that are left in any one sample */
    unsigned long most_bytes_left = 
        (cs_stats.n_loci_total - cs_stats.n_loci_read)
        * max_bytes_per_locus / cs_stats.n_threads;

    return MAX(most_bytes_left, cs_stats.n_min_absolute_bytes);
}
