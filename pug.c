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
            "-l INT      maximum length of a pileup line in bytes.  {Nuisance parameter} [100000]\n"
            "-s INT      target leaf size of a file chunk in which no more index building is done [1e7]\n"
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

int less_contig_index(const void *pa, const void *pb)
{
    const struct contig_index
        *a = (struct contig_index *)pa,
        *b = (struct contig_index *)pb;

    return strcmp(a->contig, b->contig);
}

struct contig_index *contig_buf, *contig_ptr;
unsigned num_contigs = 0;

// std::map<char *, size_t, ltstr> contig_order;
// std::map<char *, size_t, ltstr>::iterator contig_iter;

/* look up a contig and assign its index if available.  assumes c_iter
   is valid to begin with, and maintains this invariant. */
#define INIT_CONTIG(contig, c_index)                                   \
    do {                                                               \
        if (strcmp(contig_ptr->contig, contig))                        \
        {                                                              \
            struct contig_index key;                                          \
            strcpy(key.contig, contig);                                \
            contig_ptr =                                               \
                (struct contig_index *)bsearch(&key, contig_buf,              \
                                        num_contigs,                   \
                                        sizeof(struct contig_index),          \
                                        less_contig_index);            \
            if (! contig_ptr)                                          \
            {                                                          \
                fprintf(stderr,                                        \
                        "Error, Contig %s not listed "                 \
                        "in contig_order.rdb file\n", contig);         \
                exit(1);                                               \
            }                                                          \
        }                                                              \
        (c_index) = contig_ptr->index;                                 \
    } while (0)


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
        (ccmp = a->contig - b->contig) != 0
        ? ccmp
        : (a->pos - b->pos);
        
}

struct locus_range {
    struct locus_pos beg, end;
};


size_t target_leaf_size = 1e7;

/* allows mapping file offsets to loci */
struct off_index {
    struct locus_range span;
    size_t start_offset, end_offset;
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


/* find (or create) a minimal index node that contains cur, using a
   two-phase search.  the expansion phase traverses up the index tree
   until ix contains cur.  the contraction phase creates new,
   successively smaller nodes as needed, until a minimal index node
   that contains cur exists or is found.  if the root index node (with
   parent == NULL) still doesn't contain cur, return 0.  return 1 on
   success.
*/
int update_index_node(struct off_index **ixp, struct locus_pos cur, FILE *pileup_fh)
{

    /* expansion phase */
    while (! contains(*ixp, cur))
    {
        if (! (*ixp)->parent) return 1;
        *ixp = (*ixp)->parent;
    }

    /* contraction phase.  assume ix contains cur. */
    char contig[50];
    unsigned pos, ci;
    struct locus_pos midpoint_loc;
    size_t midpoint_off;
    int line_start;

    struct off_index *ix = *ixp;
    while (ix->end_offset - ix->start_offset > target_leaf_size)
    {
        /* find the midpoint */
        if (ix->left == NULL && ix->right == NULL)
        {
            fseek(pileup_fh, (ix->end_offset + ix->start_offset) / 2, SEEK_SET);
            size_t partial_pos = (size_t)ftell(pileup_fh);
            fscanf(pileup_fh, "%*[^\n]\n%n%s\t%u\t", &line_start, contig, &pos);
            INIT_CONTIG(contig, ci);
            midpoint_loc.contig = ci;
            midpoint_loc.pos = pos;
            midpoint_off = partial_pos + line_start;
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
        }

        /* traverse to appropriate child node */
        *ixp = cmp < 0 ? ix->left : ix->right;
        ix = *ixp;
    }
    return 0;
}


void free_index(struct off_index *root)
{
    if (root->left) free_index(root->left);
    if (root->right) free_index(root->right);
    free(root);
}

/* ix is the index node that contains cur.  

   1. reads pileup_fh region marked by ix.  2. finds the sub-region
   between cur and query->end. 3. writes out this range.  4. updates
   cur to the end of processed loci, which is either ix->span.end or
   query->end, whichever is less. */
void process_chunk(struct off_index *ix,
                   char *chunk_buf,
                   FILE *pileup_fh,
                   struct locus_range query,
                   struct locus_pos *cur)
{
    fseek(pileup_fh, (long)ix->start_offset, SEEK_SET);
    size_t nbytes_wanted = ix->end_offset - ix->start_offset;
    size_t nbytes_read = fread(chunk_buf, 1, nbytes_wanted, pileup_fh);
    assert(nbytes_read == nbytes_wanted);
    char *start = chunk_buf, *end, *chunk_end = chunk_buf + nbytes_read, contig[50];
    unsigned ci, pos;

    assert(chunk_buf[nbytes_read - 1] == '\n');

    /* scan forwards, updating locus and start, until */
    struct locus_pos locus = ix->span.beg;
    
    while (less_locus_pos(&locus, cur) < 0)
    {
        start = strchr(start, '\n') + 1;
        if (start == chunk_end)
        {
            *cur = ix->span.end;
            return;
        }
        int nparsed = sscanf(start, "%s\t%u\t", contig, &pos);
        assert(nparsed == 2);
        
        INIT_CONTIG(contig, ci);
        locus.contig = ci;
        locus.pos = pos;
    }

    /* there exists at least one existing locus that is not less than
       c. continue scanning forwards to find the first locus not less
       than query.end */
    end = start;
    while (less_locus_pos(&locus, &query.end) < 0)
    {
        end = strchr(end, '\n') + 1;
        if (end == chunk_end)
        {
            locus = ix->span.end;
            break;
        }
        
        int nparsed = sscanf(end, "%s\t%u\t", contig, &pos);
        assert(nparsed == 2);
        INIT_CONTIG(contig, ci);
        locus.contig = ci;
        locus.pos = pos;
    }
    /* just write to stdout */
    write(1, start, end - start);
    fflush(stdout);
    *cur = locus;

}
                   


int main_pug(int argc, char ** argv)
{
    char c;
    size_t max_pileup_line_size = 1e6;

    while ((c = getopt(argc, argv, "l:s:")) >= 0)
    {
        switch(c)
        {
        case 'l': max_pileup_line_size = (size_t)atof(optarg); break;
        case 's': target_leaf_size = (size_t)atof(optarg); break;
        default: return pug_usage(); break;
        }
    }
    if (argc - optind != 3)
        return pug_usage();

    const char *pileup_file = argv[optind];
    const char *locus_file = argv[optind + 1];
    const char *contig_order_file = argv[optind + 2];

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
    
    qsort(contig_buf, num_contigs, sizeof(contig_buf[0]), less_contig_index);

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
    fseek(pileup_fh, -max_pileup_line_size, SEEK_END);
    root->end_offset = (size_t)ftell(pileup_fh);
    root->left = root->right = root->parent = NULL;

    size_t len = fread(target_line, 1, max_pileup_line_size, pileup_fh);
    assert(len > 1);
    char *last_line = (char *)memrchr(target_line, '\n', len - 1) + 1;
    sscanf(last_line, "%s\t%u\t", contig, &root->span.end.pos);
    ++root->span.end.pos;
    INIT_CONTIG(contig, root->span.end.contig);

    fseek(pileup_fh, 0, SEEK_SET);

    /* main loop */
    q = queries;
    qend = queries + num_queries;

    struct off_index *ix = root;
    struct locus_pos cur = q->beg;

    char *chunk_buf = (char *)malloc(target_leaf_size + 1);

    while (q != qend)
    {
        update_index_node(&ix, cur, pileup_fh);
        process_chunk(ix, chunk_buf, pileup_fh, *q, &cur);

        if (less_locus_pos(&cur, &q->end) >= 0)
        {
            ++q;
            if (q != qend)
                cur = q->beg;
        }
        
    }

    free(contig_buf);
    free(chunk_buf);
    free(target_line);
    free(queries);
    free_index(root);

    fclose(pileup_fh);

    return 0;
}
