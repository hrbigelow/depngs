#define _GNU_SOURCE
#define _FILE_OFFSET_BITS 64

#include <stdio.h>
#include <unistd.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>

#include "dict.h"
#include "cache.h"
#include "pileup_bsearch.h"


int pug_usage()
{
    fprintf(stderr,
            "\nUsage: dep pug [options] sample.pileup ranges.rdb contig_order.rdb\n"
            "Options:\n\n"
            "-m INT      maximum memory to use for write buffer [1e8]\n"
            "-h FLAG     if present, interpret ranges as [begin, end) (half-inclusive)\n"
            "            if absent, interpret ranges as [begin, end] (full-inclusive) [absent]\n"
            "\n"
            "ranges.rdb has lines like:\n"
            "chr1<tab>100<tab>200\n"
            "chr1<tab>300<tab>400\n"
            "...\n"
            "It need not be in order.  It is okay if loci do not appear in sample.pileup.\n"
            "\n"
            "contig_order.rdb has lines like:\n"
            "chr1<tab>1\n"
            "chr2<tab>2\n"
            "...\n"
            "It must be complete and consistent with order of contigs in sample.pileup.\n"
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
    int half_open_ranges = 0;

    while ((c = getopt(argc, argv, "m:h")) >= 0)
    {
        switch(c)
        {
        case 'm': max_chunk_size = (size_t)atof(optarg); break;
        case 'h': half_open_ranges = 1; break;
        default: return pug_usage(); break;
        }
    }
    if (argc - optind != 3)
        return pug_usage();

    const char *pileup_file = argv[optind];
    const char *locus_file = argv[optind + 1];
    const char *contig_order_file = argv[optind + 2];


    /* the size of a chunk of file for which it is more efficient to
       read this whole chunk rather than continue to scan through the
       file */
    size_t scan_thresh_size = 1e6;
    file_bsearch_init(init_locus, scan_thresh_size);
    
    chunk_buffer = (char *)malloc(max_chunk_size);
    if (! chunk_buffer)
    {
        fprintf(stderr, "Error: Couldn't allocate %Zu bytes for write buffer\n", max_chunk_size);
        exit(1);
    }

    // 0. parse contig_order file
    FILE *contig_order_fh = fopen(contig_order_file, "r");
    if (! contig_order_fh)
    {
        fprintf(stderr, "Couldn't open contig order file %s\n", contig_order_file);
        exit(1);
    }

    char contig[1024];
    unsigned index;
    while (! feof(contig_order_fh))
    {
        fscanf(contig_order_fh, "%s\t%u\n", contig, &index);
        dict_add_item(contig, index);
    }
    fclose(contig_order_fh);
    
    dict_build();

    FILE *locus_fh = fopen(locus_file, "r");
    if (! locus_fh)
    {
        fprintf(stderr, "Couldn't open locus file %s\n", locus_file);
        exit(1);
    }

    /* 1. parse all query ranges into 'queries' and sort them */
    unsigned num_queries = 0, num_alloc = 10;
    struct locus_range *queries = 
        (struct locus_range *)malloc(num_alloc * sizeof(struct locus_range)),
        *qend,
        *q;
    
    char reformat_buf[1000];
    unsigned beg_pos, end_pos;
    while (! feof(locus_fh))
    {
        if (fscanf(locus_fh, "%s\t%u\t%u\n", contig, &beg_pos, &end_pos) != 3)
        {
            fprintf(stderr, "Error: ranges file %s has invalid format\n", locus_file);
            exit(1);
        }

        sprintf(reformat_buf, "%s\t%u\t", contig, beg_pos);
        queries[num_queries].beg = init_locus(reformat_buf);

        sprintf(reformat_buf, "%s\t%u\t", contig, end_pos);
        queries[num_queries].end = init_locus(reformat_buf);

        ++num_queries;
        ALLOC_GROW(queries, num_queries + 1, num_alloc);
    }   
    fclose(locus_fh);

    qsort(queries, num_queries, sizeof(queries[0]), less_locus_range);
    
    /* 2. edit queries to eliminate interval overlap */
    struct locus_range *p = NULL;
    for (q = queries; q != queries + num_queries - 1; ++q)
    {
        if (p && less_file_bsearch_ord(&p->end, &q->beg) > 0)
        {
            /* must be on same contig if they are overlapping and
               sorted */
            assert(p->end.hi == q->beg.hi);
            q->beg.lo = p->end.lo;
            if (q->end.lo < q->beg.lo)
                q->end.lo = q->beg.lo;
        }
        p = q;
    }


    FILE *pileup_fh = fopen(pileup_file, "r");
    if (! pileup_fh)
    {
        fprintf(stderr, "Error: Couldn't find pileup file %s\n", pileup_file);
        exit(1);
    }

    /* 3. create and initialize a root index node representing the
       entire pileup file */
    struct file_bsearch_index *root = find_root_index(pileup_fh);

    /* main loop */
    q = queries;
    qend = queries + num_queries;
    
    struct file_bsearch_index *lbeg = root, *lend;
    off_t tbeg, tend;
    while (q != qend)
    {
        lbeg = find_loose_index(lbeg, q->beg, pileup_fh);
        tbeg = off_lower_bound(lbeg, q->beg);
        lend = find_loose_index(lbeg, q->end, pileup_fh); /* lbeg first argument intentional */
        tend = half_open_ranges ? off_lower_bound(lend, q->end) : off_upper_bound(lend, q->end);

        fprintf(stderr, "Processed %zu: %zu-%zu\n", q->beg.hi, q->beg.lo, q->end.lo);
        fseeko(pileup_fh, tbeg, SEEK_SET);
        size_t bytes_to_write = tend - tbeg;
        while (bytes_to_write)
        {
            size_t chunk_bytes = MIN(bytes_to_write, max_chunk_size);
            size_t nbytes_read = fread(chunk_buffer, 1, chunk_bytes, pileup_fh);
            assert(nbytes_read == chunk_bytes);
            write(1, chunk_buffer, nbytes_read);
            bytes_to_write -= chunk_bytes;
        }
        ++q;
    }

    free(queries);
    file_bsearch_index_free(root);
    dict_free();
    file_bsearch_free();

    fclose(pileup_fh);

    return 0;
}
