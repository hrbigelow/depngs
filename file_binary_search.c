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

    struct file_bsearch_index ix = { fh, root, root };
    return ix;
}


/* find a loose-fitting index that contains cur, starting at ix, using
   only binary search.  Does not guarantee to find the tightest
   possible index because the bisection is done based on file offsets,
   not loci positions. returns NULL if the root doesn't even contain
   cur */
struct file_bsearch_node *
find_loose_index(struct file_bsearch_node *ix, struct pair_ordering cur, FILE *fh)
{
    /* expansion phase */
    while (! contains(ix, &cur))
    {
        if (! ix->parent) 
        {
            fprintf(stderr, "Error at %s: index does not contain query\n", __func__);
            exit(1);
        }
        ix = ix->parent;
    }
    /* contraction phase.  assume ix contains cur. */
    struct pair_ordering midpoint_ord;
    off_t midpoint_off, read_off = -1;
    size_t span;
    ssize_t nchars_read;
    char *line_start;
    while (1)
    {
        /* find the midpoint */
        if (read_off == -1
            && (span = ix->end_offset - ix->start_offset) <= mem_scan_threshold)
        {
            read_off = ix->start_offset;
            fseeko(fh, ix->start_offset, SEEK_SET);
            fread(mem_scan_buf, 1, span, fh);
        }
        if (ix->left == NULL && ix->right == NULL)
        {
            midpoint_off = (ix->end_offset + ix->start_offset) / 2;
            if (read_off != -1)
            {
                /* initialize midpoint_off, midpoint_ord from memory */
                line_start = strchr(mem_scan_buf + (midpoint_off - read_off), '\n') + 1;
                midpoint_off = read_off + (line_start - mem_scan_buf);
                if (midpoint_off == ix->end_offset)
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
                if (midpoint_off == ix->end_offset)
                    break;
                
                nchars_read = getline(&line_buf, &line_len, fh);
                assert(nchars_read != -1);

                line_start = line_buf;
                midpoint_ord = get_line_ord(line_start);
            }

        }
        else
        {
            midpoint_ord = ix->left ? ix->left->span_end : ix->right->span_beg;
            midpoint_off = ix->left ? ix->left->end_offset : ix->right->start_offset;
        }

        /* create a new child node as necessary */
        int cmp = cmp_pair_ordering(&cur, &midpoint_ord);
        if (cmp < 0 && ! ix->left)
        {
            ix->left = (struct file_bsearch_node *)malloc(sizeof(struct file_bsearch_node));
            struct file_bsearch_node *p = ix->left;
            p->span_beg = ix->span_beg;
            p->span_end = midpoint_ord;
            p->start_offset = ix->start_offset;
            p->end_offset = midpoint_off;
            p->parent = ix;
            p->left = p->right = NULL;
            p->span_contents = NULL;
            assert(p->start_offset < p->end_offset);
        }            
        if (cmp >= 0 && ! ix->right)
        {
            ix->right = (struct file_bsearch_node *)malloc(sizeof(struct file_bsearch_node));
            struct file_bsearch_node *p = ix->right;
            p->span_beg = midpoint_ord;
            p->span_end = ix->span_end;
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
    ix->span_contents = strndup(mem_scan_buf + (ix->start_offset - read_off),
                                ix->end_offset - ix->start_offset);
    return ix;

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
    ix->cur_node = find_loose_index(ix->cur_node, query, ix->fh);
    return off_bound_aux(ix->cur_node, query, 0);
}

/* return the smallest offset such that all lines spanning this
   query end at or before this offset. */
off_t off_upper_bound(struct file_bsearch_index *ix,
                      struct pair_ordering query)
{
    ix->cur_node = find_loose_index(ix->cur_node, query, ix->fh);
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
void file_bsearch_node_free(struct file_bsearch_node *root)
{
    if (root->left) file_bsearch_node_free(root->left);
    if (root->right) file_bsearch_node_free(root->right);
    if (root->span_contents) 
        free(root->span_contents);
    free(root);
}


void file_bsearch_index_free(struct file_bsearch_index ix)
{
    file_bsearch_node_free(ix.root);
}
