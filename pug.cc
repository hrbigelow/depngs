#include <cstdio>
#include <unistd.h>
#include <cstdlib>
#include <map>
#include <algorithm>

#include <sys/stat.h>

#include "locus_comp.h"
#include "tools.h"
#include "../samutil/file_utils.h"

int pug_usage()
{
    fprintf(stderr,
            "\nUsage: dep pug [options] sample.pileup loci_to_retrieve.rdb contig_order.rdb\n"
            "Options:\n\n"
            "-m INT      number bytes of memory to use [1e9]\n"
            "-l INT      maximum length of a pileup line in bytes.  {Nuisance parameter} [100000]\n"
            "-b INT      size of output buffer [8e6]\n"
            "\n"
            "loci_to_retrieve.rdb has lines like:\n"
            "chr1<tab>19583984\n"
            "chr1<tab>19598940\n"
            "...\n"
            "It need not be in order.  Also, it is okay if loci do not appear in the sample.pileup file\n"
            "\n"
            "contig_order.rdb has lines of <contig><tab><ordering>\n"
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


std::map<char const*, size_t, ltstr> contig_order;
std::map<char const*, size_t, ltstr>::iterator contig_iter;

/* look up a contig and assign its index if available.  assumes c_iter
   is valid to begin with, and maintains this invariant. */
#define INIT_CONTIG(contig, c_index)                                   \
    do {                                                               \
        if (strcmp((*contig_iter).first, contig))                      \
        {                                                              \
            contig_iter = contig_order.find(contig);                   \
            if (contig_iter == contig_order.end())                     \
            {                                                          \
                fprintf(stderr,                                        \
                        "Error, Contig %s not listed "                 \
                        "in contig_order.rdb file\n", contig);         \
                exit(1);                                               \
            }                                                          \
        }                                                              \
        (c_index) = (*contig_iter).second;                             \
    } while (0)


struct locus_pos { 
    unsigned contig, pos;
};


int less_locus_pos(const void *pa, const void *pb)
{
    const struct locus_pos
        *a = (struct locus_pos *)pa,
        *b = (struct locus_pos *)pb;

    return 
        (a->contig - b->contig) || (a->pos - b->pos);
        
}

struct locus_range {
    locus_pos beg, end;
};


size_t target_leaf_size = 1e7;

/* allows mapping file offsets to loci */
struct off_index {
    locus_range span;
    size_t start_offset, end_offset;
    off_index *left, *right, *parent;
};


int less_locus_range(const void *pa, const void *pb)
{
    const struct locus_range
        *a = (struct locus_range *)pa,
        *b = (struct locus_range *)pb;

    return 
        less_locus_pos(&a->beg, &b->beg) || less_locus_pos(&a->end, &b->end);
}


int contains(struct off_index *ix, locus_pos loc)
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
int update_index_node(struct off_index **ixp, locus_pos cur, FILE *pileup_fh)
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
            fseek(pileup_fh, (long)(ix->end_offset - ix->start_offset), SEEK_SET);
            size_t partial_pos = (size_t)ftell(pileup_fh);
            fscanf(pileup_fh, "%*[^\n]\n%n%s\t%u\t", &line_start, contig, &pos);
            INIT_CONTIG(contig, ci);
            midpoint_loc = { ci, pos };
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
            p->span = { ix->span.beg, midpoint_loc };
            p->start_offset = ix->start_offset;
            p->end_offset = midpoint_off;
            p->parent = ix;
            p->left = p->right = NULL;
        }            
        else if (cmp >= 0 && ! ix->right)
        {
            ix->right = (struct off_index *)malloc(sizeof(struct off_index));
            struct off_index *p = ix->right;
            p->span = { midpoint_loc, ix->span.end };
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
   query->end, whichever is less.  if  */
void process_chunk(struct off_index *ix,
                   struct locus_pos *cur,
                   char *chunk_buf,
                   FILE *pileup_fh,
                   struct locus_range **query)
{
    fseek(pileup_fh, (long)ix->start_offset, SEEK_SET);
    size_t nbytes_wanted = ix->end_offset - ix->start_offset;
    size_t nbytes_read = fread(chunk_buf, 1, nbytes_wanted, pileup_fh);
    assert(nbytes_read == nbytes_wanted);
    char *start = chunk_buf, *end = start + nbytes_read, contig[50];
    unsigned ci, pos;
    /* scan forwards */
    struct locus_pos locus = ix->span.beg;
    while (less_locus_pos(&locus, cur) < 0)
    {
        int nparsed = sscanf(start, "%s\t%u\t", contig, &pos);
        assert(nparsed == 2);
        start = strchr(start, '\n') + 1;

        INIT_CONTIG(contig, ci);
        locus = { ci, pos };
    }

    /* scan backwards in the chunk to find greatest line less than
       query->end */
    locus = ix->span.end;
    while (less_locus_pos(&locus, &(*query)->end) >= 0)
    {
        end = (char *)memrchr(chunk_buf, '\n', end - chunk_buf);
        int nparsed = sscanf(end, "%s\t%u\t", contig, &pos);
        assert(nparsed == 2);
        INIT_CONTIG(contig, ci);
        locus = { ci, pos };
    }
    /* just write to stdout */
    write(1, start, end - start);
    *cur = locus;

    if (less_locus_pos(cur, (*query)->end) == 0)
        /* we have output all of the loci in this query */
        ++(*query);
}
                   


int main_pug(int argc, char ** argv)
{
    char c;
    size_t max_mem = 1024l * 1024l * 1024l;
    size_t max_pileup_line_size = 1e6;
    size_t outbuf_size = 8e6;

    while ((c = getopt(argc, argv, "m:l:b:s:")) >= 0)
    {
        switch(c)
        {
        case 'm': max_mem = (size_t)atof(optarg); break;
        case 'l': max_pileup_line_size = (size_t)atof(optarg); break;
        case 'b': outbuf_size = (size_t)atof(optarg); break;
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
    FILE * contig_order_fh = open_if_present(contig_order_file, "r");
    while (! feof(contig_order_fh))
    {
        char *contig = (char *)malloc(100);
        size_t order;
        fscanf(contig_order_fh, "%s\t%zu\n", contig, & order);
        contig_order.insert(std::make_pair(contig, order));
    }
    fclose(contig_order_fh);
    assert(! contig_order.empty());
    
    /* invariant: contig_iter is always valid */
    contig_iter = contig_order.begin();
    
    FILE *locus_fh = open_if_present(locus_file, "r");

    /* 1. parse all query ranges into 'queries' and sort them */
    char contig[200];
    unsigned num_queries = 0, num_alloc = 10;
    struct locus_range *queries = 
        (struct locus_range *)malloc(num_alloc * sizeof(struct locus_range)),
        *qend,
        *q;
    
    while (fscanf(locus_fh, "%s\t%u\t%u\n", contig, &q->beg.pos, &q->end.pos) == 3)
    {
        INIT_CONTIG(contig, q->beg.contig);
        q->end.contig = q->beg.contig; /* query ranges are assumed to be from the same contig */
        ++q;
        ++num_queries;
        ALLOC_GROW(q, num_queries, num_alloc);
    }   
    fclose(locus_fh);

    qend = queries + num_queries;

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
    
    FILE *pileup_fh = open_if_present(pileup_file, "r");

    char *target_line = (char *)malloc(max_pileup_line_size + 1);

    /* 3. create and initialize a root index node representing the entire pileup file */
    struct off_index *root = (struct off_index *)malloc(sizeof(struct off_index));
    char contig[50];
    root->start_offset = 0;
    fscanf("%s\t%u\t", contig, &root->beg->pos, pileup_fh);
    INIT_CONTIG(contig, root->beg->contig);
    fseek(pileup_fh, -max_pileup_line_size, SEEK_END);
    root->end_offset = (size_t)ftell(pileup_fh);
    root->left = root->right = root->parent = NULL;

    size_t len = fread(target_line, 1, max_pileup_line_size, pileup_fh);
    assert(len > 1);
    char *last_line = memrchr(target_line, '\n', len - 1) + 1;
    sscanf("%s\t%u\t", contig, &root->end->pos, last_line);
    ++root->end->pos;
    INIT_CONTIG(contig, root->end->contig);

    fseek(pileup_fh, 0, SEEK_SET);

    /* main loop */
    struct off_index *ix = root;
    struct locus_pos cur = q.beg;

    char *chunk_buf = (char *)malloc(target_leaf_size + 1);

    while (q != qend)
    {
        update_index_node(ix, &cur, pileup_fh);
        process_chunk(&ix, &cur, chunk_buf, pileup_fh, &q);
    }

    for (contig_iter = contig_order.begin();
         contig_iter != contig_order.end(); 
         ++contig_iter)
        free((*contig_iter).first);

    free(chunk_buf);
    free(target_line);
    free(queries);
    free_index(root);

    fclose(pileup_fh);

    return 0;
}
