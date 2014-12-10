/* functions for doing binary search on a file. */
#include "file_binary_search.h"

#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <stdint.h>

#define CMP(a, b) ((a) < (b) ? -1 : (a) > (b) ? 1 : 0)
#define MIN(a,b) ((a) < (b) ? (a) : (b))
#define MAX(a,b) ((a) < (b) ? (b) : (a))

/* defines the ordering logic for the file_bsearch_ord structure. */
int less_file_bsearch_ord(const void *pa, const void *pb)
{
    const struct file_bsearch_ord
        *a = (struct file_bsearch_ord *)pa,
        *b = (struct file_bsearch_ord *)pb;

    int acmp;
    return
        (acmp = CMP(a->hi, b->hi)) != 0
        ? acmp
        : CMP(a->lo, b->lo);
}


/* returns 1 if ix contains ord, 0 otherwise */
int contains(struct file_bsearch_index *ix, struct file_bsearch_ord *ord)
{
    return (less_file_bsearch_ord(&ix->span.beg, ord) <= 0
            && less_file_bsearch_ord(ord, &ix->span.end) < 0);
}

/* parses an input line to get its ordinal structure */
static struct file_bsearch_ord (*get_line_ord)(const char *line);

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
void file_bsearch_init(struct file_bsearch_ord (*_get_line_ord)(const char *line),
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

#if 0
/* generate a root index that spans the whole file */
struct file_bsearch_index *find_root_index(FILE *fh)
{
    struct file_bsearch_index *root = 
        (struct file_bsearch_index *) malloc(sizeof(struct file_bsearch_index));

    if (! get_line_ord)
    {
        fprintf(stderr, "file_binary_search: error, you didn't call file_besearch_init()\n");
        exit(1);
    }

    root->start_offset = 0;
    fseeko(fh, 0, SEEK_SET);
    getline(&line_buf, &line_len, fh);
    root->span.beg = get_line_ord(line_buf);

    /* scan backwards from the end until we find a newline */
    fseeko(fh, 0, SEEK_END);
    root->end_offset = ftello(fh);

    /* we assume that the file ends in a newline. this is the position
       of the penultimate newline. */
    char *penul_newline = NULL;
    off_t cur = root->end_offset;
    while (! penul_newline && cur != 0)
    {
        cur -= MIN(line_len, cur);
        fseeko(fh, cur, SEEK_SET);
        size_t remain = MIN(line_len, root->end_offset - cur);
        (void)fread(line_buf, 1, remain, fh);
        penul_newline = (char *)memrchr(line_buf, '\n', remain - 1);
    }
    root->span.end = get_line_ord(penul_newline + 1);
    root->span.end.lo++; /* this is the next position of ord after the
                            penultimate. */

    root->left = root->right = root->parent = NULL;
    root->span_contents = NULL;
    return root;
}
#endif


static struct file_bsearch_ord min_ord = { 0, 0 };
static struct file_bsearch_ord max_ord = { SIZE_MAX, SIZE_MAX };

/* generate a root index that spans the whole file */
struct file_bsearch_index *find_root_index(FILE *fh)
{
    struct file_bsearch_index *root = 
        (struct file_bsearch_index *) malloc(sizeof(struct file_bsearch_index));

    if (! get_line_ord)
    {
        fprintf(stderr, "file_binary_search: error, you didn't call file_besearch_init()\n");
        exit(1);
    }
    root->span.beg = min_ord;
    root->span.end = max_ord;
    root->start_offset = 0;
    fseeko(fh, 0, SEEK_END);
    root->end_offset = ftello(fh);
    root->left = root->right = root->parent = NULL;
    root->span_contents = NULL;

    return root;
}


/* find a loose-fitting index that contains cur, starting at ix, using
   only binary search.  Does not guarantee to find the tightest
   possible index because the bisection is done based on file offsets,
   not loci positions. returns NULL if the root doesn't even contain
   cur */
struct file_bsearch_index *
find_loose_index(struct file_bsearch_index *ix, struct file_bsearch_ord cur, FILE *fh)
{
    /* expansion phase */
    while (! contains(ix, &cur))
    {
        if (! ix->parent) return NULL;
        ix = ix->parent;
    }
    /* contraction phase.  assume ix contains cur. */
    struct file_bsearch_ord midpoint_ord;
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
            midpoint_ord = ix->left ? ix->left->span.end : ix->right->span.beg;
            midpoint_off = ix->left ? ix->left->end_offset : ix->right->start_offset;
        }

        /* create a new child node as necessary */
        int cmp = less_file_bsearch_ord(&cur, &midpoint_ord);
        if (cmp < 0 && ! ix->left)
        {
            ix->left = (struct file_bsearch_index *)malloc(sizeof(struct file_bsearch_index));
            struct file_bsearch_index *p = ix->left;
            p->span.beg = ix->span.beg;
            p->span.end = midpoint_ord;
            p->start_offset = ix->start_offset;
            p->end_offset = midpoint_off;
            p->parent = ix;
            p->left = p->right = NULL;
            p->span_contents = NULL;
            assert(p->start_offset < p->end_offset);
        }            
        if (cmp >= 0 && ! ix->right)
        {
            ix->right = (struct file_bsearch_index *)malloc(sizeof(struct file_bsearch_index));
            struct file_bsearch_index *p = ix->right;
            p->span.beg = midpoint_ord;
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
    ix->span_contents = strndup(mem_scan_buf + (ix->start_offset - read_off),
                                ix->end_offset - ix->start_offset);
    return ix;

}

inline
off_t off_bound_aux(const struct file_bsearch_index *ix, struct file_bsearch_ord query, int cmp)
{
    struct file_bsearch_ord cur = ix->span.beg;
    char
        *start = ix->span_contents,
        *end = start + (ix->end_offset - ix->start_offset);

    while (less_file_bsearch_ord(&cur, &query) < cmp)
    {
        start = strchr(start, '\n') + 1;
        if (start == end)
        {
            cur = ix->span.end;
            break;
        }
        cur = get_line_ord(start);
    }
    return ix->start_offset + (start - ix->span_contents);
}


off_t off_lower_bound(const struct file_bsearch_index *ix, struct file_bsearch_ord query)
{
    return off_bound_aux(ix, query, 0);
}

off_t off_upper_bound(const struct file_bsearch_index *ix, struct file_bsearch_ord query)
{
    return off_bound_aux(ix, query, 1);
}


/* free the index tree */
void file_bsearch_index_free(struct file_bsearch_index *root)
{
    if (root->left) file_bsearch_index_free(root->left);
    if (root->right) file_bsearch_index_free(root->right);
    if (root->span_contents) 
        free(root->span_contents);
    free(root);
}

