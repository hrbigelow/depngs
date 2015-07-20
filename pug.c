#define _GNU_SOURCE
#define _FILE_OFFSET_BITS 64

#include <stdio.h>
#include <unistd.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>

#include "cache.h"
#include "locus.h"
#include "file_binary_search.h"

int pug_usage()
{
    fprintf(stderr,
            "\nUsage: dep pug [options] sample.pileup ranges.rdb index.fai\n"
            "Options:\n\n"
            "-m INT      maximum memory to use for write buffer [1e8]\n"
            "\n"
            "ranges.rdb has lines like:\n"
            "chr1<tab>100<tab>200\n"
            "chr1<tab>300<tab>400\n"
            "...\n"
            "It need not be in order.  Each contig in ranges.rdb must appear in\n"
            "index.fai file.  Ranges are interpreted as half-open: 100-200 means\n"
            "loci in [100, 200), that is, up to but not including 200\n"
            "\n"
            "index.fai is the fasta index file produced from samtools faidx\n"
            "\n"
            );
    return 1;
}

#define MIN(a,b) ((a) < (b) ? (a) : (b))

static size_t max_chunk_size = 1e8;
static char *chunk_buffer;

int main_pug(int argc, char ** argv)
{
    char c;
    while ((c = getopt(argc, argv, "m")) >= 0) {
        switch(c) {
        case 'm': max_chunk_size = (size_t)atof(optarg); break;
        default: return pug_usage(); break;
        }
    }
    if (argc - optind != 3)
        return pug_usage();

    const char *pileup_file = argv[optind];
    const char *locus_file = argv[optind + 1];
    const char *faidx_file = argv[optind + 2];

    locus_global_init(faidx_file);

    /* the size of a chunk of file for which it is more efficient to
       read this whole chunk rather than continue to scan through the
       file */
    size_t scan_thresh_size = 1e6;
    file_bsearch_init(init_locus, scan_thresh_size);
    
    chunk_buffer = (char *)malloc(max_chunk_size);
    if (! chunk_buffer) {
        fprintf(stderr, "Error: Couldn't allocate %Zu bytes for write buffer\n", max_chunk_size);
        exit(1);
    }

    /* parse fasta index file */
    faidx_t *fai = fai_load(faidx_file);
    if (fai == NULL) {
        fprintf(stderr, "%s:%i: Error: couldn't parse or load fasta index file %s\n",
                __FILE__, __LINE__, faidx_file);
        exit(1);
    }

    unsigned n_queries;
    unsigned long num_total_loci;
    struct pair_ordering_range *queries =
        parse_query_ranges(locus_file, fai, &n_queries, &num_total_loci);

    struct pair_ordering_range *q = queries, *qend = q + n_queries;

    /* main loop */
    struct file_bsearch_index ix = file_bsearch_make_index(pileup_file);
    struct pair_ordering cur_beg, cur_end;
    while (q != qend) {
        fprintf(stderr, "Processed %zu: %zu-%zu\n", q->beg.hi, q->beg.lo, q->end.lo);
        size_t bytes_to_write = range_to_size(&ix, q->beg, q->end);

        cur_beg = q->beg;

        while (bytes_to_write) {
            cur_end = bytes_to_write <= max_chunk_size
                ? q->end
                : size_to_range(&ix, cur_beg, max_chunk_size);
            size_t nbytes_read = read_range(&ix, cur_beg, cur_end, chunk_buffer);
            (void)write(1, chunk_buffer, nbytes_read);
            bytes_to_write -= nbytes_read;
            cur_beg = cur_end;
        }
        ++q;
        if (ix.n_nodes > 1000)
            (void)file_bsearch_range_free(&ix, 
                                          (struct pair_ordering){ 0, 0 },
                                          (struct pair_ordering)
                                          { ix.root->span_end.hi, ix.root->span_end.lo - 1 });
    }
    free(queries);
    file_bsearch_index_free(ix);
    locus_global_free();

    return 0;
}
