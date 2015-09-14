#include "chunk_strategy.h"

#include <stdlib.h>

struct chunk_strategy cs_stats;

#define MAX(a,b) ((a) < (b) ? (b) : (a))

/* in order to simplify the interface, we choose a reasonable default
   value for this. it may be off by a factor of 10 or so, but that
   won't matter for large input. */
#define DEFAULT_BYTES_PER_LOCUS 10

void
chunk_strategy_init(unsigned n_files, unsigned n_threads,
                    const char *locus_range_file,
                    const char *fasta_file,
                    unsigned long bytes_zone2,
                    unsigned long bytes_zone3)
{
    cs_stats.n_files = n_files;
    cs_stats.n_threads = n_threads;
    cs_stats.bytes_zone2 = bytes_zone2;
    cs_stats.bytes_zone3 = bytes_zone3;
    cs_stats.query_regions =
        parse_locus_ranges(locus_range_file,
                           fasta_file,
                           &cs_stats.n_query_regions,
                           &cs_stats.n_total_loci);

    if (cs_stats.n_query_regions > 1000)
    cs_stats.n_all_bytes_read = calloc(n_files, sizeof(cs_stats.n_all_bytes_read[0]));
    cs_stats.n_all_loci_read = 0;

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
    cs_stats.total_span = span;
    const struct contig_region *qlo, *qhi;
    cs_stats.n_total_loci =
        find_intersecting_span(cs_stats.query_regions,
                               cs_stats.query_regions + cs_stats.n_query_regions,
                               span, &qlo, &qhi);
    cs_stats.cur_pos = span.beg;
}


/* return the estimated maximum number of bytes that we should attempt
   to parse per thread.  divide each zone into N_PIECES_PER_THREAD *
   n_threads.  This is enough pieces so that each thread is exposed to
   10 pieces on average, allowing the stochastic differences in piece
   size to average out.  */
#define N_PIECES_PER_THREAD 10
unsigned long
cs_max_bytes_wanted()
{
    /* calculate how many loci are left */
    struct contig_span subset = { cs_stats.cur_pos, cs_stats.total_span.end };
    const struct contig_region *qlo, *qhi;
    unsigned long n_loci_left = 
        find_intersecting_span(cs_stats.query_regions,
                               cs_stats.query_regions + cs_stats.n_query_regions,
                               subset,
                               &qlo, &qhi);
    
    /* estimate the minimum number of bytes_left in any sample. if
       there are lots of query regions, there will be inefficient
       parsing because BAM scanning has a finest resolution of one
       tile of 16384 bp. So, we assume we will be parsing an entire
       BGZF block (16000 bases on average) for one locus. */
    unsigned long max_bytes_per_locus = 0;
    if (cs_stats.n_all_loci_read == 0)
        max_bytes_per_locus =
            cs_stats.n_query_regions > 100
            ? 16000
            : 10;
    else {
        unsigned long max_bytes_read = 0;
        unsigned f;
        for (f = 0; f != cs_stats.n_files; ++f)
            max_bytes_read = MAX(max_bytes_read, cs_stats.n_all_bytes_read[f]);
        max_bytes_per_locus = max_bytes_read / cs_stats.n_all_loci_read;
    }

    /* the most bytes that are left in any one sample */
    unsigned long most_bytes_left = n_loci_left * max_bytes_per_locus;
    unsigned long n_pieces = cs_stats.n_threads * N_PIECES_PER_THREAD;

    /* return the size of the zone, divided by the number of threads */
    if (most_bytes_left > cs_stats.bytes_zone2 + cs_stats.bytes_zone3) {
        /* in zone 1 */
        unsigned long n_est_bytes_zone1 = cs_stats.n_total_loci * max_bytes_per_locus;
        return n_est_bytes_zone1 / n_pieces;
    } else if (most_bytes_left > cs_stats.bytes_zone3) {
        /* in zone 2 */
        return cs_stats.bytes_zone2 / n_pieces;
    } else {
        /* in zone 3 */
        return cs_stats.bytes_zone3 / n_pieces;
    }
}
