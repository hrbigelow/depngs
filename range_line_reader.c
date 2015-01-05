#define _GNU_SOURCE
#define _FILE_OFFSET_BITS 64

#include "range_line_reader.h"

#include <stddef.h>
#include <stdlib.h>
#include <assert.h>
#include <string.h>
#include <stdio.h>

#define MAX(a, b) ((a) > (b) ? (a) : (b))
#define MIN(a, b) ((a) < (b) ? (a) : (b))

#define MIN_PAIR_ORD(a, b) (cmp_pair_ordering(&(a), &(b)) < 0 ? (a) : (b))

/* Read the same amount of the specified ranges in multiple files in
   tandem. par specifies the number of buffers */
void range_line_reader(void *par, struct managed_buf *bufs)
{
    struct range_line_reader_par *rr = par;
    size_t n_ix = rr->n_ix;

    off_t max_span;
    off_t *span = malloc(n_ix * sizeof(off_t));

    char
        **write_ptr = malloc(n_ix * sizeof(char *)),
        **write_end = malloc(n_ix * sizeof(char *));

    size_t i;
    for (i = 0; i != n_ix; ++i)
    {
        write_ptr[i] = bufs[i].buf;
        write_end[i] = bufs[i].buf + bufs[i].alloc;
    }
    ptrdiff_t space_left = write_end[0] - write_ptr[0], space_tmp;
    struct pair_ordering tmp_pos, trunc_pos = max_pair_ord;

    /* must re-set on each call */
    rr->new_query = 1;
    
    while (rr->new_query && rr->q != rr->qend)
    {
        max_span = 0;
        for (i = 0; i != n_ix; ++i)
        {
            span[i] = range_to_size(&rr->ix[i], rr->q->beg, rr->q->end);
            max_span = MAX(max_span, span[i]);
        }

        /* fill the buffer as much as possible with the next query range.
           afterwards, q now points to the next range to retrieve. */
        if (max_span <= space_left)
        {
            for (i = 0; i != n_ix; ++i)
            {
                write_ptr[i] += read_range(&rr->ix[i], rr->q->beg, rr->q->end, 
                                           write_ptr[i]);
                space_tmp = write_end[i] - write_ptr[i];
                space_left = MIN(space_left, space_tmp);
            }
            ++rr->q;
            rr->new_query = 1;
        }
        else
        {
            for (i = 0; i != n_ix; ++i)
            {
                tmp_pos = size_to_range(&rr->ix[i], rr->q->beg, space_left);
                trunc_pos = MIN_PAIR_ORD(trunc_pos, tmp_pos);
            }
            /* this assumes all buffers are equal size */
            if (space_left == bufs[0].alloc
                && cmp_pair_ordering(&rr->q->beg, &trunc_pos) == 0)
            {
                fprintf(stderr, 
                        "%s: Couldn't fit the line at %Zu:%Zu "
                        "into the allotted space of %Zu\n",
                        __func__, rr->q->beg.hi, rr->q->beg.lo, space_left);
                exit(1);
            }
            
            for (i = 0; i != n_ix; ++i)
            {
                write_ptr[i] += read_range(&rr->ix[i], rr->q->beg, trunc_pos, 
                                           write_ptr[i]);
                space_tmp = write_end[i] - write_ptr[i];
                space_left = MIN(space_left, space_tmp);
            }
            rr->q->beg = trunc_pos;
            rr->new_query = 0;
        }
    }
    for (i = 0; i != n_ix; ++i)
        bufs[i].size = write_ptr[i] - bufs[i].buf;

    /* free resources */
    for (i = 0; i != n_ix; ++i)
        (void)file_bsearch_node_range_free(rr->ix[i].root,
                                           min_pair_ord,
                                           rr->q == rr->qend ? max_pair_ord : rr->q->beg);

    free(span);
    free(write_ptr);
    free(write_end);
}
