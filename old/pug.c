#define _GNU_SOURCE
#define _FILE_OFFSET_BITS 64

#include <stdio.h>
#include <unistd.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>

#include "cache.h"
#include "locus.h"
#include "locus_range.h"
#include "file_binary_search.h"
#include "fasta.h"

int pug_usage()
{
    fprintf(stderr,
            "\nUsage: dep pug [options] sample.pileup ranges.rdb ref.fa.fai\n"
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
            "ref.fa.fai is the fasta index file.  it is only used to define the number\n"
            "and order of contigs.  it must match that used to produce the pileup file\n"
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
    const char *fai_file = argv[optind + 2];

    /* strangely, we need to first strip off the .fai extension */
    char fasta_file[1000];
    unsigned fl = strlen(fai_file);
    if (strcmp(fai_file + fl - 4, ".fai") != 0) {
        fprintf(stderr, "%s:%u Error: fasta index file %s must have the '.fai' extension\n",
                __FILE__, __LINE__, fai_file);
        exit(1);
    }
    strncpy(fasta_file, fai_file, fl - 4);
    fasta_file[fl] = '\0';
    locus_init(fasta_file);

    /* the size of a chunk of file for which it is more efficient to
       read this whole chunk rather than continue to scan through the
       file */
    size_t scan_thresh_size = 1e6;
    file_bsearch_init(parse_pileup_locus, scan_thresh_size);
    
    chunk_buffer = (char *)malloc(max_chunk_size);
    if (! chunk_buffer) {
        fprintf(stderr, "Error: Couldn't allocate %Zu bytes for write buffer\n", max_chunk_size);
        exit(1);
    }

    unsigned n_queries;
    unsigned long num_total_loci;
    struct contig_region *queries =
        parse_locus_ranges(locus_file, fasta_file, &n_queries, &num_total_loci);

    struct contig_region *q = queries, *qend = q + n_queries;

    /* main loop.  annoyingly, we are working with contig_region (tid,
       beg, end) and pair_ordering (hi, lo). tid is like 'hi', and lo
       is like beg or end) */
    struct file_bsearch_index ix = file_bsearch_make_index(pileup_file);
    struct pair_ordering cur_beg, cur_end;
    while (q != qend) {
        fprintf(stderr, "Processed %s: %u-%u\n", 
                fasta_get_contig(q->tid), q->beg, q->end);
        struct pair_ordering 
            pbeg = { q->tid, q->beg },
            pend = { q->tid, q->end };
            size_t bytes_to_write = range_to_size(&ix, pbeg, pend);

        cur_beg = pbeg;

        while (bytes_to_write) {
            cur_end = bytes_to_write <= max_chunk_size
                ? pend
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
    locus_free();
    free(queries);
    file_bsearch_index_free(ix);

    return 0;
}
