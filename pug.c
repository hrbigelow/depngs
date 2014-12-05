#define _GNU_SOURCE

#include <stdio.h>
#include <unistd.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>

#include "dict.h"
#include "cache.h"

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

/* GOTCHA!  Do NOT subtract two unsigned's and expect it to convert to
   a signed quantity */
// #define CMP(a, b) ((int)(a) - (int)(b))
// #define CMP(a, b) ((a) - (b))


size_t max_chunk_size = 1e8;
char *chunk_buffer;

char *pileup_read_buf;
char *pileup_line_buf;

off_t pileup_read_off; /* kept in synch with pileup_read_buf.  -1
                          means pileup_read_buf is undefined.  other
                          values mean that pileup_read_buf holds the
                          contents of pileup_fh at that offset. */

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
void init_locus(const char *line, 
                struct locus_pos *loc)
{
    char contig[50];
    unsigned pos;
    int nparsed = sscanf(line, "%s\t%u\t", contig, &pos);
    assert(nparsed == 2);
    long ix;
    if ((ix = dict_search(contig)) >= 0)
        loc->contig = (unsigned)ix;
    else
        MISSING_CONTIG(contig);

    loc->pos = pos;

}

    
int less_locus_pos(const void *pa, const void *pb)
{
    const struct locus_pos
        *a = (struct locus_pos *)pa,
        *b = (struct locus_pos *)pb;

    int ccmp;
    return
        (ccmp = CMP(a->contig, b->contig)) != 0
        ? ccmp
        : CMP(a->pos, b->pos);
        
}

struct locus_range {
    struct locus_pos beg, end;
};


/* allows mapping file offsets to loci */
struct off_index {
    struct locus_range span;
    off_t start_offset, end_offset;
    struct off_index *left, *right, *parent;
    char *span_contents; /* if non-null, will contain contents of
                            file */
};


int less_locus_range(const void *pa, const void *pb)
{
    const struct locus_range
        *a = (struct locus_range *)pa,
        *b = (struct locus_range *)pb;

    int bcmp;
    return 
        (bcmp = less_locus_pos(&a->beg, &b->beg)) != 0
        ? bcmp
        : less_locus_pos(&a->end, &b->end);
}


int contains(struct off_index *ix, struct locus_pos loc)
{
    return (less_locus_pos(&ix->span.beg, &loc) <= 0
            && less_locus_pos(&loc, &ix->span.end) < 0);
}


/* find a loose-fitting index that contains cur, starting at ix, using
   only binary search.  Does not guarantee to find the tightest
   possible index because the bisection is done based on file offsets,
   not loci positions. returns NULL if the root doesn't even contain
   cur */
struct off_index *
find_loose_index(struct off_index *ix, struct locus_pos cur, FILE *pileup_fh)
{
    /* expansion phase */
    while (! contains(ix, cur))
    {
        if (! ix->parent) return NULL;
        ix = ix->parent;
    }
    /* contraction phase.  assume ix contains cur. */
    char contig[50];
    unsigned pos;
    struct locus_pos midpoint_loc;
    off_t midpoint_off;
    int rval;
    size_t dummy, span;
    ssize_t nchars_read;
    pileup_read_off = -1;
    char *line_start;
    while (1)
    {
        /* find the midpoint */
        if (pileup_read_off == -1
            && (span = ix->end_offset - ix->start_offset) <= max_chunk_size)
        {
            pileup_read_off = ix->start_offset;
            fseeko(pileup_fh, ix->start_offset, SEEK_SET);
            fread(pileup_read_buf, 1, span, pileup_fh);
            /* buf_end = pileup_read_buf + span; */
        }
        if (ix->left == NULL && ix->right == NULL)
        {
            midpoint_off = (ix->end_offset + ix->start_offset) / 2;
            if (pileup_read_off != -1)
            {
                /* initialize midpoint_off, midpoint_loc from memory */
                line_start = strchr(pileup_read_buf + (midpoint_off - pileup_read_off), '\n') + 1;
                midpoint_off = pileup_read_off + (line_start - pileup_read_buf);
                if (midpoint_off == ix->end_offset)
                    break;
                rval = sscanf(line_start, "%s\t%u\t", contig, &pos);
                assert(rval == 2);
            }
            else
            {
                /* initialize midpoint_off, midpoint_loc from file */
                fseeko(pileup_fh, midpoint_off, SEEK_SET);
                nchars_read = getline(&pileup_line_buf, &dummy, pileup_fh);
                assert(nchars_read != -1);
                
                midpoint_off += nchars_read;
                if (midpoint_off == ix->end_offset)
                    break;
                
                rval = fscanf(pileup_fh, "%s\t%u\t", contig, &pos);
                assert(rval == 2);
            }

            long cix = dict_search(contig);
            if (cix < 0)
                MISSING_CONTIG(contig);
            midpoint_loc.contig = (unsigned)cix;
            midpoint_loc.pos = pos;
            
        }
        else
        {
            midpoint_loc = ix->left ? ix->left->span.end : ix->right->span.beg;
            midpoint_off = ix->left ? ix->left->end_offset : ix->right->start_offset;
        }

        /* create a new child node as necessary */
        int cmp = less_locus_pos(&cur, &midpoint_loc);
        if (cmp < 0 && ! ix->left)
        {
            ix->left = (struct off_index *)malloc(sizeof(struct off_index));
            struct off_index *p = ix->left;
            p->span.beg = ix->span.beg;
            p->span.end = midpoint_loc;
            p->start_offset = ix->start_offset;
            p->end_offset = midpoint_off;
            p->parent = ix;
            p->left = p->right = NULL;
            p->span_contents = NULL;
            assert(p->start_offset < p->end_offset);
        }            
        if (cmp >= 0 && ! ix->right)
        {
            ix->right = (struct off_index *)malloc(sizeof(struct off_index));
            struct off_index *p = ix->right;
            p->span.beg = midpoint_loc;
            p->span.end = ix->span.end;
            p->start_offset = midpoint_off;
            p->end_offset = ix->end_offset;
            p->parent = ix;
            p->left = p->right = NULL;
            p->span_contents = NULL;
            assert(p->start_offset < p->end_offset);
        }

        /* traverse to appropriate child node */
        ix = cmp < 0 ? ix->left : ix->right;
    }
    ix->span_contents = strndup(pileup_read_buf + (ix->start_offset - pileup_read_off),
                                ix->end_offset - ix->start_offset);
    return ix;

}


/* Scan forward locus by locus, creating a tightest index containing
   cur. */
#define OFF_BOUND(WHICH, OP)                                            \
    off_t off_ ## WHICH ## _bound(const struct off_index *ix,           \
                                  struct locus_pos cur)                 \
    {                                                                   \
        struct locus_pos locus = ix->span.beg;                          \
        char                                                            \
            *start = ix->span_contents,                                 \
            *end = start + (ix->end_offset - ix->start_offset);         \
                                                                        \
        while (less_locus_pos(&locus, &cur) OP 0)                       \
        {                                                               \
            start = strchr(start, '\n') + 1;                            \
            if (start == end)                                           \
            {                                                           \
                locus = ix->span.end;                                   \
                break;                                                  \
            }                                                           \
            init_locus(start, &locus);                                   \
        }                                                               \
        return ix->start_offset + (start - ix->span_contents);          \
    }                                                                   \
    

/* searches ix for the file offset of the lower bound of cur */
OFF_BOUND(lower, <)

/* searches ix for the file offset of the upper bound of cur */
OFF_BOUND(upper, <=)

void free_index(struct off_index *root)
{
    if (root->left) free_index(root->left);
    if (root->right) free_index(root->right);
    if (root->span_contents) 
        free(root->span_contents);
    free(root);
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
    
    while (fscanf(locus_fh, "%s\t%u\t%u\n", contig, 
                  &queries[num_queries].beg.pos, 
                  &queries[num_queries].end.pos) == 3)
    {
        cix = dict_search(contig);
        if (cix < 0)
            MISSING_CONTIG(contig);

        queries[num_queries].beg.contig = (unsigned)cix;
        queries[num_queries].end.contig = (unsigned)cix;
        ++num_queries;
        ALLOC_GROW(queries, num_queries + 1, num_alloc);
    }   
    fclose(locus_fh);

    qsort(queries, num_queries, sizeof(queries[0]), less_locus_range);
    
    /* 2. edit queries to eliminate interval overlap */
    struct locus_range *p = NULL;
    for (q = queries; q != queries + num_queries - 1; ++q)
    {
        if (p && less_locus_pos(&p->end, &q->beg) > 0)
        {
            /* must be on same contig if they are overlapping and
               sorted */
            assert(p->end.contig == q->beg.contig);
            q->beg.pos = p->end.pos;
            if (q->end.pos < q->beg.pos)
                q->end.pos = q->beg.pos;
        }
        p = q;
    }

    pileup_read_buf = (char *)malloc(max_chunk_size);
    if (! pileup_read_buf)
    {
        fprintf(stderr, "Error: Couldn't allocate a read buffer of size %zu\n", max_chunk_size);
        exit(1);
    }

    pileup_line_buf = (char *)malloc(max_pileup_line_size);
    if (! pileup_line_buf)
    {
        fprintf(stderr, "Error: Couldn't allocate a line buffer of size %zu\n", 
                max_pileup_line_size);
        exit(1);
    }


    FILE *pileup_fh = fopen(pileup_file, "r");
    if (! pileup_fh)
    {
        fprintf(stderr, "Error: Couldn't find pileup file %s\n", pileup_file);
        exit(1);
    }

    // setbuf(pileup_fh, NULL);

    char *target_line = (char *)malloc(max_pileup_line_size + 1);

    /* 3. create and initialize a root index node representing the entire pileup file */
    struct off_index *root = (struct off_index *) malloc(sizeof(struct off_index));

    root->start_offset = 0;
    fscanf(pileup_fh, "%s\t%u\t", contig, &root->span.beg.pos);
    cix = dict_search(contig);
    if (cix < 0)
        MISSING_CONTIG(contig);
    root->span.beg.contig = (unsigned)cix;

    fseeko(pileup_fh, (off_t)-max_pileup_line_size, SEEK_END);
    size_t len = fread(target_line, 1, max_pileup_line_size, pileup_fh);
    root->end_offset = ftello(pileup_fh);
    assert(len > 1);
    
    root->left = root->right = root->parent = NULL;
    root->span_contents = NULL;

    char *last_line = (char *)memrchr(target_line, '\n', len - 1) + 1;
    init_locus(last_line, &root->span.end);
    root->span.end.pos++;

    free(target_line);

    fseeko(pileup_fh, 0, SEEK_SET);

    /* main loop */
    q = queries;
    qend = queries + num_queries;

    struct off_index *lbeg = root, *lend;
    off_t tbeg, tend;
    while (q != qend)
    {
        lbeg = find_loose_index(lbeg, q->beg, pileup_fh);
        tbeg = off_lower_bound(lbeg, q->beg);
        lend = find_loose_index(lbeg, q->end, pileup_fh); /* lbeg first argument intentional */
        tend = off_upper_bound(lend, q->end);

        fprintf(stderr, "Processed %u: %u-%u\n", q->beg.contig, q->beg.pos, q->end.pos);
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
    free_index(root);
    dict_free();

    fclose(pileup_fh);

    return 0;
}
