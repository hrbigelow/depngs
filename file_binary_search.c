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

/* size of buffer to hold a section of the file that is deemed small enough to
   be fully loaded rather than continue traversing with fseeko */
static size_t mem_scan_threshold;


/* load buf with enough content from fh, and return the location in
   buf the nearest line start of a complete line to the left of the
   current position in the file.  buf is allocated with an initial
   'sz'. */
char *rfind_linestart(char *buf, size_t buf_sz, FILE *fh)
{
    char *nl = NULL;
    off_t max_sz = ftello(fh);
    off_t off = MIN(buf_sz, max_sz);
    do
    {
        buf = realloc(buf, off);
        fseeko(fh, -off, SEEK_CUR);
        fread(buf, 1, off, fh);
        nl = memrchr(buf, '\n', off);
        if (nl) nl = memrchr(buf, '\n', nl - buf);
        off *= 2;
        off = MIN(off, max_sz);
    }
    while (! nl && off != max_sz);
    return nl ? nl + 1 : buf;
}


/* initialize resources */
void file_bsearch_init(get_line_ord_t _get_line_ord,
                       size_t _mem_scan_threshold)
{
    get_line_ord = _get_line_ord;
    mem_scan_threshold = _mem_scan_threshold;
}

/* generate a root index that spans the whole file */
struct file_bsearch_index file_bsearch_make_index(const char *file)
{
    struct file_bsearch_node *root = 
        malloc(sizeof(struct file_bsearch_node));

    FILE *fh = fopen(file, "r");
    if (! fh)
    {
        fprintf(stderr, "%s: couldn't open file %s\n", __func__, file);
        exit(1);
    }
    if (! get_line_ord)
    {
        fprintf(stderr, "file_binary_search: error, you didn't call file_besearch_init()\n");
        exit(1);
    }
    setvbuf(fh, NULL, _IONBF, 0);

    struct file_bsearch_index ix = { 
        .fh = fh,
        .mem_scan_buf = malloc(mem_scan_threshold),
        .line_buf = malloc(1e5),
        .line_len = 1e5,
        .root = root, 
        .cur_node = root, 
        .n_nodes = 1
    };

    int nchars_read = getline(&ix.line_buf, &ix.line_len, ix.fh);
    assert(nchars_read != -1);

    root->span_beg = get_line_ord(ix.line_buf);
    root->start_offset = 0;

    fseeko(ix.fh, 0, SEEK_END);
    root->end_offset = ftello(fh);
    unsigned scan_sz = MIN(root->end_offset - root->start_offset, 1e3);
    char *end_buf = malloc(scan_sz);
    
    char *line_start = rfind_linestart(end_buf, scan_sz, ix.fh);
    root->span_end = get_line_ord(line_start);
    root->span_end.lo++;
    free(end_buf);

    root->left = root->right = root->parent = NULL;
    root->span_contents = NULL;

    return ix;
}


/* find a loose-fitting index that contains cur, or the root node if
   not found.  starting at ix->cur_node, using only binary search.
   Does not guarantee to find the tightest possible index because the
   bisection is done based on file offsets, not loci
   positions. updates ix->cur_node to point to this new node.  */
void
find_loose_index(struct file_bsearch_index *ix, struct pair_ordering cur, FILE *fh)
{
    /* expansion phase */
    struct file_bsearch_node *nd = ix->cur_node;
    while (! contains(nd, &cur) && nd->parent)
    {
        /* if (! nd->parent)  */
        /* { */
        /*     fprintf(stderr, "Error at %s: index does not contain query\n", __func__); */
        /*     exit(1); */
        /* } */
        nd = nd->parent;
    }

    /* contraction phase.  assume nd contains cur. */
    struct pair_ordering midpoint_ord;

    /* read_off = -1 indicates that mem_scan_buf has not been
       initialized.  if it remains -1 after the first test, then it
       we go into file initialization mode. */
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
            (void)fread(ix->mem_scan_buf, 1, span, fh);
        }
        if (nd->left == NULL && nd->right == NULL)
        {
            midpoint_off = (nd->end_offset + nd->start_offset) / 2;
            if (read_off != -1)
            {
                /* initialize midpoint_off, midpoint_ord from memory */
                line_start = strchr(ix->mem_scan_buf + (midpoint_off - read_off), '\n') + 1;
                midpoint_off = read_off + (line_start - ix->mem_scan_buf);
                if (midpoint_off == nd->end_offset)
                    break;
                
                midpoint_ord = get_line_ord(line_start);
            }
            else
            {
                /* initialize midpoint_off, midpoint_ord from file */
                fseeko(fh, midpoint_off, SEEK_SET);
                nchars_read = getline(&ix->line_buf, &ix->line_len, fh);
                assert(nchars_read != -1);
                
                midpoint_off += nchars_read;
                if (midpoint_off == nd->end_offset)
                    break;
                
                nchars_read = getline(&ix->line_buf, &ix->line_len, fh);
                assert(nchars_read != -1);

                line_start = ix->line_buf;
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
            assert(nd->start_offset < midpoint_off);
            ix->n_nodes++;
            nd->left = malloc(sizeof(struct file_bsearch_node));
            *nd->left = (struct file_bsearch_node){
                nd->span_beg, midpoint_ord,
                nd->start_offset, midpoint_off,
                NULL, NULL, nd, NULL
            };
        }            
        if (cmp >= 0 && ! nd->right)
        {
            assert(midpoint_off < nd->end_offset);
            ix->n_nodes++;
            nd->right = malloc(sizeof(struct file_bsearch_node));
            *nd->right = (struct file_bsearch_node){
                midpoint_ord, nd->span_end,
                midpoint_off, nd->end_offset,
                NULL, NULL, nd, NULL
            };
        }

        /* traverse to appropriate child node */
        nd = cmp < 0 ? nd->left : nd->right;
    }
    if (read_off != -1)
    {
        /* mem_scan_buf has been initialized.  safe to populate span_contents */
        nd->span_contents = 
            strndup(ix->mem_scan_buf + (nd->start_offset - read_off), span);
    }
    else
    {
        nd->span_contents = malloc(span);
        fseeko(ix->fh, nd->start_offset, SEEK_SET);
        fread(nd->span_contents, 1, span, ix->fh);
    }
    ix->cur_node = nd;
}


/* search forward in nd->span_contents for the offset which is the
   start (cmp == 0) or end (cmp == 1) of first line greater than
   query. if query < nd->span_beg, return nd->start_offset.  if query
   >= nd->span_end, return nd->end_offset */
off_t off_bound_aux(struct file_bsearch_node *nd,
                    struct pair_ordering query,
                    int cmp)
{
    struct pair_ordering cur = nd->span_beg;

    char
        *start = nd->span_contents,
        *end = start + (nd->end_offset - nd->start_offset);

    while (cmp_pair_ordering(&cur, &query) < cmp)
    {
        start = strchr(start, '\n') + 1;
        if (start == end)
        {
            cur = nd->span_end;
            break;
        }
        cur = get_line_ord(start);
    }
    return nd->start_offset + (start - nd->span_contents);
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
    off_t off_cur = off_lower_bound(ix, beg) + (off_t)size;

    if (off_cur > ix->root->end_offset)
        return (struct pair_ordering){ SIZE_MAX, SIZE_MAX };

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
    free(ix.mem_scan_buf);
    if (ix.line_buf) free(ix.line_buf);
    if (ix.fh) fclose(ix.fh);
    (void)file_bsearch_node_free(ix.root);
}

/* free all nodes that are contained in [beg, end)*/
size_t node_range_free(struct file_bsearch_node *node,
                       struct pair_ordering beg,
                       struct pair_ordering end)
{
    struct pair_ordering mid;
    int cmp_beg, cmp_end;
    size_t n_freed = 0;

    if (! node) return 0;
    if ((cmp_beg = cmp_pair_ordering(&beg, &node->span_beg)) <= 0
        && (cmp_end = cmp_pair_ordering(&node->span_end, &end)) <= 0)
    {
        struct file_bsearch_node *parent = node->parent;
        n_freed = file_bsearch_node_free(node);
        if (parent && parent->left == node) parent->left = NULL;
        if (parent && parent->right == node) parent->right = NULL;
        return n_freed;
    }
    mid = node->left ? node->left->span_end 
        : node->right ? node->right->span_beg
        : beg;
    
    if (cmp_pair_ordering(&beg, &mid) < 0)
        n_freed += node_range_free(node->left, beg, end);

    if (cmp_pair_ordering(&mid, &end) < 0)
        n_freed += node_range_free(node->right, beg, end);

    return n_freed;
}


/* frees all nodes in [beg, end) range.  sets cur_node to root. */
size_t file_bsearch_range_free(struct file_bsearch_index *ix,
                               struct pair_ordering beg,
                               struct pair_ordering end)
{
    ix->cur_node = ix->root;
    return node_range_free(ix->root, beg, end);
}
