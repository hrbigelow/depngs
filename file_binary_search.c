/* functions for doing binary search on a file. */
#include "file_binary_search.h"

#include <stdlib.h>
#include <string.h>
#include <assert.h>

#define MIN(a,b) ((a) < (b) ? (a) : (b))
#define MAX(a,b) ((a) < (b) ? (b) : (a))


/* returns 1 if ix contains ord, 0 otherwise */
int contains(struct file_bsearch_node *ix, struct pair_ordering *ord)
{
    return (cmp_pair_ordering(&ix->span_beg, ord) <= 0
            && cmp_pair_ordering(ord, &ix->span_end) < 0);
}

/* parses an input line to get its ordinal structure */
static get_line_ord_t get_line_ord;

/* buffer to hold a section of the file that is deemed small enough to
   be fully loaded rather than continue traversing with fseeko */
static char *mem_scan_buf;
static size_t mem_scan_threshold;


/* a very conservative estimate for maximum line size that we expect.
   This will grow as needed if longer lines are encountered. */
static size_t line_len = 1e5; 

/* reallocated as needed if longer lines are encountered. */
static char *line_buf; 

/* initialize resources */
void file_bsearch_init(get_line_ord_t _get_line_ord,
                       size_t _mem_scan_threshold)
{
    get_line_ord = _get_line_ord;
    mem_scan_threshold = _mem_scan_threshold;
    mem_scan_buf = (char *)malloc(mem_scan_threshold);
    line_buf = (char *)malloc(line_len);
}

/* release resources */
void file_bsearch_free()
{
    if (mem_scan_buf) 
    {
        free(mem_scan_buf);
        mem_scan_buf = NULL;
        mem_scan_threshold = 0;
    }
    if (line_buf) 
    {
        free(line_buf);
        line_buf = NULL;
        line_len = 0;
    }
    get_line_ord = NULL;
}


/* generate a root index that spans the whole file */
struct file_bsearch_index file_bsearch_make_index(FILE *fh)
{
    struct file_bsearch_node *root = 
        (struct file_bsearch_node *) malloc(sizeof(struct file_bsearch_node));

    if (! get_line_ord)
    {
        fprintf(stderr, "file_binary_search: error, you didn't call file_besearch_init()\n");
        exit(1);
    }
    root->span_beg = min_pair_ord;
    root->span_end = max_pair_ord;
    root->start_offset = 0;
    fseeko(fh, 0, SEEK_END);
    root->end_offset = ftello(fh);
    root->left = root->right = root->parent = NULL;
    root->span_contents = NULL;

    struct file_bsearch_index ix = { 
        .fh = fh, 
        .root = root, 
        .cur_node = root, 
        .n_nodes = 1, 
    };

    return ix;
}


/* find a loose-fitting index that contains cur, starting at
   ix->cur_node, using only binary search.  Does not guarantee to find
   the tightest possible index because the bisection is done based on
   file offsets, not loci positions. updates ix->cur_node to point to
   this new node. updates */
void
find_loose_index(struct file_bsearch_index *ix, struct pair_ordering cur, FILE *fh)
{
    /* expansion phase */
    struct file_bsearch_node *nd = ix->cur_node;
    while (! contains(nd, &cur))
    {
        if (! nd->parent) 
        {
            fprintf(stderr, "Error at %s: index does not contain query\n", __func__);
            exit(1);
        }
        nd = nd->parent;
    }

    /* contraction phase.  assume nd contains cur. */
    struct pair_ordering midpoint_ord;
    off_t midpoint_off, read_off = -1;
    size_t span;
    ssize_t nchars_read;
    char *line_start;
    while (1)
    {
        /* find the midpoint */
        if (read_off == -1
            && (span = nd->end_offset - nd->start_offset) <= mem_scan_threshold)
        {
            read_off = nd->start_offset;
            fseeko(fh, nd->start_offset, SEEK_SET);
            fread(mem_scan_buf, 1, span, fh);
        }
        if (nd->left == NULL && nd->right == NULL)
        {
            midpoint_off = (nd->end_offset + nd->start_offset) / 2;
            if (read_off != -1)
            {
                /* initialize midpoint_off, midpoint_ord from memory */
                line_start = strchr(mem_scan_buf + (midpoint_off - read_off), '\n') + 1;
                midpoint_off = read_off + (line_start - mem_scan_buf);
                if (midpoint_off == nd->end_offset)
                    break;
                
                midpoint_ord = get_line_ord(line_start);
            }
            else
            {
                /* initialize midpoint_off, midpoint_ord from file */
                fseeko(fh, midpoint_off, SEEK_SET);
                nchars_read = getline(&line_buf, &line_len, fh);
                assert(nchars_read != -1);
                
                midpoint_off += nchars_read;
                if (midpoint_off == nd->end_offset)
                    break;
                
                nchars_read = getline(&line_buf, &line_len, fh);
                assert(nchars_read != -1);

                line_start = line_buf;
                midpoint_ord = get_line_ord(line_start);
            }

        }
        else
        {
            midpoint_ord = nd->left ? nd->left->span_end : nd->right->span_beg;
            midpoint_off = nd->left ? nd->left->end_offset : nd->right->start_offset;
        }

        /* create a new child node as necessary */
        int cmp = cmp_pair_ordering(&cur, &midpoint_ord);
        if (cmp < 0 && ! nd->left)
        {
            nd->left = (struct file_bsearch_node *)malloc(sizeof(struct file_bsearch_node));
            ix->n_nodes++;
            struct file_bsearch_node *p = nd->left;
            p->span_beg = nd->span_beg;
            p->span_end = midpoint_ord;
            p->start_offset = nd->start_offset;
            p->end_offset = midpoint_off;
            p->parent = nd;
            p->left = p->right = NULL;
            p->span_contents = NULL;
            assert(p->start_offset < p->end_offset);
        }            
        if (cmp >= 0 && ! nd->right)
        {
            nd->right = (struct file_bsearch_node *)malloc(sizeof(struct file_bsearch_node));
            ix->n_nodes++;
            struct file_bsearch_node *p = nd->right;
            p->span_beg = midpoint_ord;
            p->span_end = nd->span_end;
            p->start_offset = midpoint_off;
            p->end_offset = nd->end_offset;
            p->parent = nd;
            p->left = p->right = NULL;
            p->span_contents = NULL;
            assert(p->start_offset < p->end_offset);
        }

        /* traverse to appropriate child node */
        nd = cmp < 0 ? nd->left : nd->right;
    }
    nd->span_contents = strndup(mem_scan_buf + (nd->start_offset - read_off),
                                nd->end_offset - nd->start_offset);
    assert(nd->end_offset - nd->start_offset < 10000);
    ix->cur_node = nd;
}

inline
off_t off_bound_aux(const struct file_bsearch_node *ix,
                    struct pair_ordering query,
                    int cmp)
{
    struct pair_ordering cur = ix->span_beg;
    char
        *start = ix->span_contents,
        *end = start + (ix->end_offset - ix->start_offset);

    while (cmp_pair_ordering(&cur, &query) < cmp)
    {
        start = strchr(start, '\n') + 1;
        if (start == end)
        {
            cur = ix->span_end;
            break;
        }
        cur = get_line_ord(start);
    }
    return ix->start_offset + (start - ix->span_contents);
}

/* together, these two functions give you the smallest chunk of file
   containing all lines that span query */

/* return the largest offset such that all lines spanning query start
   at or after this offset. */
off_t off_lower_bound(struct file_bsearch_index *ix,
                      struct pair_ordering query)
{
    find_loose_index(ix, query, ix->fh);
    return off_bound_aux(ix->cur_node, query, 0);
}

/* return the smallest offset such that all lines spanning this
   query end at or before this offset. */
off_t off_upper_bound(struct file_bsearch_index *ix,
                      struct pair_ordering query)
{
    find_loose_index(ix, query, ix->fh);
    return off_bound_aux(ix->cur_node, query, 1);
}


/* return an estimated upper size limit to contain the logical range
   [beg, end) */
size_t range_to_size(struct file_bsearch_index *ix,
                     struct pair_ordering beg,
                     struct pair_ordering end)
{
    return off_lower_bound(ix, end) - off_lower_bound(ix, beg);
}


/* return an estimated highest upper bound logical position that uses
   <= size bytes */
struct pair_ordering size_to_range(struct file_bsearch_index *ix,
                                   struct pair_ordering beg,
                                   size_t size)
{
    off_t 
        off_cur = off_lower_bound(ix, beg) + (off_t)size;

    if (off_cur > ix->root->end_offset)
        return max_pair_ord;

    fseeko(ix->fh, off_cur, SEEK_SET);
    off_t off = MIN(1000, size); /* reasonable initial estimate for a line length */
    char *buf = malloc(off + 1), *newline = NULL;
    
    while (! newline && off != size)
    {
        buf = realloc(buf, off + 1);
        fseeko(ix->fh, -off, SEEK_CUR);
        fread(buf, 1, off, ix->fh);
        newline = memrchr(buf, '\n', off);
        if (newline) newline = memrchr(buf, '\n', newline - buf);
        off *= 2;
        off = MIN(off, size);
    }
    struct pair_ordering end = newline ? get_line_ord(newline + 1) : beg;
    free(buf);
    return end;
}


/* read a logical range into a buffer.  caller must ensure buf is
   allocated. */
size_t read_range(struct file_bsearch_index *ix,
                  struct pair_ordering beg,
                  struct pair_ordering end,
                  char *buf)
{
    off_t off_beg = off_lower_bound(ix, beg);
    off_t off_end = off_lower_bound(ix, end);
    fseeko(ix->fh, off_beg, SEEK_SET);
    size_t n_read = fread(buf, 1, off_end - off_beg, ix->fh);
    assert(n_read == off_end - off_beg);
    return n_read;
}



/* free the index tree */
size_t file_bsearch_node_free(struct file_bsearch_node *root)
{
    size_t c = 0;
    if (root->left) c += file_bsearch_node_free(root->left);
    if (root->right) c += file_bsearch_node_free(root->right);
    if (root->span_contents) 
        free(root->span_contents);
    free(root);
    c++;
    return c;
}


void file_bsearch_index_free(struct file_bsearch_index ix)
{
    (void)file_bsearch_node_free(ix.root);
}

/* free all nodes that are contained in [beg, end)*/
size_t file_bsearch_node_range_free(struct file_bsearch_node *node,
                                    struct pair_ordering beg,
                                    struct pair_ordering end)
{
    struct pair_ordering mid;
    int cmp_beg, cmp_end;
    size_t n_freed = 0;

    if (! node) return 0;
    if ((cmp_beg = cmp_pair_ordering(&beg, &node->span_beg)) <= 0
        && (cmp_end = cmp_pair_ordering(&end, &node->span_end)) <= 0)
    {
        struct file_bsearch_node *parent = node;
        n_freed = file_bsearch_node_free(node);
        if (parent && parent->left == node) parent->left = NULL;
        if (parent && parent->right == node) parent->right = NULL;
        return n_freed;
    }
    mid = node->left ? node->left->span_end 
        : node->right ? node->right->span_beg
        : beg;
    
    if (cmp_pair_ordering(&beg, &mid) < 0)
        n_freed += file_bsearch_node_range_free(node->left, beg, end);

    if (cmp_pair_ordering(&mid, &end) < 0)
        n_freed += file_bsearch_node_range_free(node->right, beg, end);

    return n_freed;
}

