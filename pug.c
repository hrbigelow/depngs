#define _GNU_SOURCE

#include <stdio.h>
#include <unistd.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>

#include "dict.h"
#include "cache.h"
#include "file_binary_search.h"

#define _FILE_OFFSET_BITS 64

int pug_usage()
{
    fprintf(stderr,
            "\nUsage: dep pug [options] sample.pileup loci_to_retrieve.rdb contig_order.rdb\n"
            "Options:\n\n"
            "-m INT      maximum memory to use for write buffer [1e8]\n"
            "-l INT      maximum length of a pileup line in bytes.  {Nuisance parameter} [100000]\n"
            "\n"
            "loci_to_retrieve.rdb has lines like:\n"
            "chr1<tab>100<tab>200\n"
            "chr1<tab>300<tab>400\n"
            "...\n"
            "This will retrieve loci from 100 to 199, and 300 to 399.\n"
            "It need not be in order.  Also, it is okay if loci do not appear in the sample.pileup file\n"
            "\n"
            "contig_order.rdb has lines of\n"
            "chr1<tab>1\n"
            "chr2<tab>2\n"
            "...\n"
            "It must be consistent and complete with the orderings of the contigs mentioned in\n"
            "all pileup input files\n"
            );
    return 1;
}

#define MIN(a,b) ((a) < (b) ? (a) : (b))
#define MAX(a,b) ((a) < (b) ? (b) : (a))

/* Do not pass arguments that are evaluated */
#define CMP(a, b) ((a) < (b) ? -1 : (a) > (b) ? 1 : 0)

size_t max_chunk_size = 1e8;
char *chunk_buffer;

#define MISSING_CONTIG(CTG)                         \
    do {                                            \
    fprintf(stderr, "Error at %s: %u. "             \
            "Contig %s not found in dictionary.\n", \
            __FILE__, __LINE__, (contig));          \
    exit(1);                                        \
} while (0)
      
      
struct locus_pos { 
    unsigned contig, pos;
};
        

/* initialize a locus from a character line */
struct file_bsearch_ord init_locus(const char *line)
{
    char contig[200];
    unsigned pos;
    struct file_bsearch_ord o;
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


struct locus_range {
    struct file_bsearch_ord beg, end;
};

          
int less_locus_range(const void *pa, const void *pb)
{
    const struct locus_range
        *a = (struct locus_range *)pa,
        *b = (struct locus_range *)pb;

    int bcmp;
    return 
        (bcmp = less_file_bsearch_ord(&a->beg, &b->beg)) != 0
        ? bcmp
        : less_file_bsearch_ord(&a->end, &b->end);
}


int main_pug(int argc, char ** argv)
{
    char c;
    size_t max_pileup_line_size = 1e6;

    while ((c = getopt(argc, argv, "l:m:")) >= 0)
    {
        switch(c)
        {
        case 'l': max_pileup_line_size = (size_t)atof(optarg); break;
        case 'm': max_chunk_size = (size_t)atof(optarg); break;
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
    long cix;
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
    while (fscanf(locus_fh, "%s\t%u\t%u\n", contig, &beg_pos, &end_pos) == 3)
    {
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
        tend = off_upper_bound(lend, q->end);

        fprintf(stderr, "Processed %zu: %zu-%zu\n", q->beg.hi, q->beg.lo, q->end.lo);
        /* fseeko(pileup_fh, tbeg.start_offset, SEEK_SET); */
        /* size_t bytes_to_write = tend.start_offset - tbeg.start_offset; */
        /* while (bytes_to_write) */
        /* { */
        /*     size_t chunk_bytes = MIN(bytes_to_write, max_chunk_size); */
        /*     size_t nbytes_read = fread(chunk_buffer, 1, chunk_bytes, pileup_fh); */
        /*     assert(nbytes_read == chunk_bytes); */
        /*     write(1, chunk_buffer, nbytes_read); */
        /*     bytes_to_write -= chunk_bytes; */
        /* } */
        ++q;
    }

    free(queries);
    file_bsearch_index_free(root);
    dict_free();
    file_bsearch_free();

    fclose(pileup_fh);

    return 0;
}
