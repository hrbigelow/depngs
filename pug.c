#define _GNU_SOURCE

#include <stdio.h>
#include <unistd.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>

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



/* taken from git source code */
#define alloc_nr(x) (((x)+16)*3/2)

/*
 * Realloc the buffer pointed at by variable 'x' so that it can hold
 * at least 'nr' entries; the number of entries currently allocated
 * is 'alloc', using the standard growing factor alloc_nr() macro.
 *
 * DO NOT USE any expression with side-effect for 'x', 'nr', or 'alloc'.
 */
#define ALLOC_GROW(x, nr, alloc)                    \
    do {                                            \
        if ((nr) > alloc) {                         \
            if (alloc_nr(alloc) < (nr))             \
                alloc = (nr);                       \
            else                                    \
                alloc = alloc_nr(alloc);            \
            x = realloc((x), alloc * sizeof(*(x))); \
        }                                           \
    } while (0)


struct contig_index {
    char contig[50];
    unsigned index;
};

int less_contig_index_c(const void *pa, const void *pb)
{
    const struct contig_index
        *a = (struct contig_index *)pa,
        *b = (struct contig_index *)pb;

    return strcmp(a->contig, b->contig);
}



struct contig_index *contig_buf, *contig_ptr;
unsigned num_contigs = 0;

size_t max_chunk_size = 1e8;
char *chunk_buffer;

/* look up a contig and assign its index if available.  assumes c_iter
   is valid to begin with, and maintains this invariant. */
#define INIT_CONTIG(contig, c_index)                                    \
    do {                                                                \
        if (strcmp(contig_ptr->contig, contig))                         \
        {                                                               \
            struct contig_index key;                                    \
            strcpy(key.contig, contig);                                 \
            contig_ptr =                                                \
                (struct contig_index *)                                 \
                bsearch(&key, contig_buf,                               \
                        num_contigs,                                    \
                        sizeof(struct contig_index),                    \
                        less_contig_index_c);                           \
            if (! contig_ptr)                                           \
            {                                                           \
                fprintf(stderr,                                         \
                        "Error, Contig %s not listed "                  \
                        "in contig_order.rdb file\n", contig);          \
                exit(1);                                                \
            }                                                           \
        }                                                               \
        (c_index) = contig_ptr->index;                                  \
    } while (0)


/* initialize a locus from a character line */
#define INIT_LOCUS(line, locus)                                   \
    do {                                                          \
        char contig[50];                                          \
        unsigned pos;                                             \
        int nparsed = sscanf(line, "%s\t%u\t", contig, &pos);     \
        assert(nparsed == 2);                                     \
        INIT_CONTIG(contig, (locus).contig);                      \
        (locus).pos = pos;                                        \
    } while (0)                                                   \

    
struct locus_pos { 
    unsigned contig, pos;
};


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
    long start_offset, end_offset;
    struct off_index *left, *right, *parent;
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
   not loci positions. */
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
    long midpoint_off = ix->end_offset;
    int rval;

    while (1)
    {
        /* find the midpoint */
        if (ix->left == NULL && ix->right == NULL)
        {
            fseek(pileup_fh, (ix->end_offset + ix->start_offset) / 2, SEEK_SET);
            int nextchar = fgetc(pileup_fh);

            if ((char)nextchar != '\n')
            {
                ungetc(nextchar, pileup_fh);
                fscanf(pileup_fh, "%*[^\n]\n"); /* throw away extra lines */
            }
            /* now, pileup_fh is queued up to the start of a line */

            midpoint_off = ftell(pileup_fh);
            if (midpoint_off == ix->end_offset)
                break;
            

            rval = fscanf(pileup_fh, "%s\t%u\t", contig, &pos);
            assert(rval == 2);
            /* assert(strncmp(contig_ptr->contig, "chr", 3) == 0); */
            INIT_CONTIG(contig, midpoint_loc.contig);
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
            assert(p->start_offset < p->end_offset);
        }

        /* traverse to appropriate child node */
        ix = cmp < 0 ? ix->left : ix->right;
        /* char c; */
        /* if (ix->start_offset != 0) */
        /* { */
        /*     fseek(pileup_fh, ix->start_offset - 1, SEEK_SET); */
        /*     fscanf(pileup_fh, "%c", &c); */
        /*     assert(c == '\n'); */
        /* } */
        /* if (ix->end_offset != 259707606333) */
        /* { */
        /*     fseek(pileup_fh, ix->end_offset - 1, SEEK_SET); */
        /*     fscanf(pileup_fh, "%c", &c); */
        /*     assert(c == '\n'); */
        /* } */
        
    }
    return ix;

}




/* Scan forward locus by locus, creating a tightest index containing
   cur.  this will be an orphan index node with no parent, and will
   not attach itself to ix */
struct off_index
get_tight_index(const struct off_index *ix, struct locus_pos cur, FILE *pileup_fh)
{
    struct off_index t;
    struct locus_pos locus = ix->span.beg;
    size_t nbytes_wanted = ix->end_offset - ix->start_offset;
    char *buf = (char *)malloc(nbytes_wanted), *start = buf, *end = buf + nbytes_wanted;
    
    
    fseek(pileup_fh, ix->start_offset, SEEK_SET);
    size_t nbytes_read = fread(buf, 1, nbytes_wanted, pileup_fh);
    assert(nbytes_read == nbytes_wanted);

    while (less_locus_pos(&locus, &cur) < 0)
    {
        start = strchr(start, '\n') + 1;
        if (start == end)
        {
            locus = ix->span.end;
            break;
        }
        INIT_LOCUS(start, locus);
    }
    t.span.beg = locus;
    t.start_offset = ix->start_offset + (start - buf);

    while (less_locus_pos(&locus, &cur) == 0)
    {
        start = strchr(start, '\n') + 1;
        if (start == end)
        {
            locus = ix->span.end;
            break;
        }
        INIT_LOCUS(start, locus);
    }
    t.span.end = locus;
    t.end_offset = ix->start_offset + (start - buf);

    free(buf);
    return t;
}

void free_index(struct off_index *root)
{
    if (root->left) free_index(root->left);
    if (root->right) free_index(root->right);
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
    unsigned num_contig_alloc = 10;
    contig_buf = (struct contig_index *)malloc(num_contig_alloc * sizeof(struct contig_index));
    contig_ptr = contig_buf;

    FILE *contig_order_fh = fopen(contig_order_file, "r");
    if (! contig_order_fh)
    {
        fprintf(stderr, "Couldn't open contig order file %s\n", contig_order_file);
        exit(1);
    }

    while (! feof(contig_order_fh))
    {
        fscanf(contig_order_fh, "%s\t%u\n", 
               contig_buf[num_contigs].contig, 
               &contig_buf[num_contigs].index);
        ++num_contigs;
        ALLOC_GROW(contig_buf, num_contigs + 1, num_contig_alloc);
    }
    fclose(contig_order_fh);
    assert(num_contigs != 0);
    
    qsort(contig_buf, num_contigs, sizeof(contig_buf[0]), less_contig_index_c);

    /* invariant: contig_iter is always valid */
    contig_ptr = contig_buf;

    FILE *locus_fh = fopen(locus_file, "r");
    if (! locus_fh)
    {
        fprintf(stderr, "Couldn't open locus file %s\n", locus_file);
        exit(1);
    }

    /* 1. parse all query ranges into 'queries' and sort them */
    char contig[50];
    unsigned num_queries = 0, num_alloc = 10;
    struct locus_range *queries = 
        (struct locus_range *)malloc(num_alloc * sizeof(struct locus_range)),
        *qend,
        *q;
    
    while (fscanf(locus_fh, "%s\t%u\t%u\n", contig, 
                  &queries[num_queries].beg.pos, 
                  &queries[num_queries].end.pos) == 3)
    {
        INIT_CONTIG(contig, queries[num_queries].beg.contig);

        /* query ranges are assumed to be from the same contig */
        queries[num_queries].end.contig = queries[num_queries].beg.contig; 
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
    
    FILE *pileup_fh = fopen(pileup_file, "r");
    if (! pileup_fh)
    {
        fprintf(stderr, "Couldn't find pileup file %s\n", pileup_file);
        exit(1);
    }

    char *target_line = (char *)malloc(max_pileup_line_size + 1);

    /* 3. create and initialize a root index node representing the entire pileup file */
    struct off_index *root = (struct off_index *) malloc(sizeof(struct off_index));

    root->start_offset = 0;
    fscanf(pileup_fh, "%s\t%u\t", contig, &root->span.beg.pos);
    INIT_CONTIG(contig, root->span.beg.contig);
    fseek(pileup_fh, -(long)max_pileup_line_size, SEEK_END);
    size_t len = fread(target_line, 1, max_pileup_line_size, pileup_fh);
    root->end_offset = ftell(pileup_fh);
    assert(len > 1);
    
    root->left = root->right = root->parent = NULL;

    char *last_line = (char *)memrchr(target_line, '\n', len - 1) + 1;
    sscanf(last_line, "%s\t%u\t", contig, &root->span.end.pos);
    ++root->span.end.pos;
    INIT_CONTIG(contig, root->span.end.contig);

    free(target_line);

    fseek(pileup_fh, 0, SEEK_SET);

    /* main loop */
    q = queries;
    qend = queries + num_queries;

    struct off_index *lbeg = root, *lend, tbeg, tend;

    while (q != qend)
    {
        lbeg = find_loose_index(lbeg, q->beg, pileup_fh);
        tbeg = get_tight_index(lbeg, q->beg, pileup_fh);
        lend = find_loose_index(lbeg, q->end, pileup_fh); /* lbeg first argument intentional */
        tend = get_tight_index(lend, q->end, pileup_fh);

        fprintf(stderr, "Processed %s: %Zu-%Zu\n", q->beg.contig, q->beg.pos, q->end.pos);
        /* fseek(pileup_fh, tbeg.start_offset, SEEK_SET); */
        /* size_t bytes_to_write = tend.start_offset - tbeg.start_offset; */
        /* while (bytes_to_write) */
        /* { */
        /*     size_t chunk_bytes = MIN(bytes_to_write, max_chunk_size); */
        /*     size_t nbytes_read = fread(chunk_buffer, 1, chunk_bytes, pileup_fh); */
        /*     assert(nbytes_read == chunk_bytes); */
        /*     write(1, chunk_buffer, nbytes_read); */
        /*     bytes_to_write -= chunk_bytes; */
        /* }         */
        ++q;
    }

    free(queries);
    free_index(root);

    fclose(pileup_fh);

    return 0;
}
