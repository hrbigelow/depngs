#define _GNU_SOURCE
#define _FILE_OFFSET_BITS 64

#include "range_line_reader.h"

#include <stddef.h>
#include <stdlib.h>
#include <assert.h>
#include <string.h>
#include <stdio.h>

#if 0
/* Conforms to the thread_queue_reader_t expectation. Manages buf.
   Assumes line-based input. */
void range_line_reader(void *par, struct managed_buf *buf)
{
    struct range_line_reader_par *rr = par;
    off_t span;
    struct pair_ordering trunc_pos;

    assert(buf->alloc > 0);
    char *write_ptr = buf->buf, *write_end = buf->buf + buf->alloc;
    ptrdiff_t space_left = write_end - write_ptr;

    /* must re-set on each call */
    rr->new_query = 1;
    
    while (rr->new_query && rr->q != rr->qend)
    {
        span = range_to_size(&rr->ix[0], rr->q->beg, rr->q->end);
        
        /* fill the buffer as much as possible with the next query range.
           afterwards, q now points to the next range to retrieve. */
        if (span <= space_left)
        {
            write_ptr += read_range(&rr->ix[0], rr->q->beg, rr->q->end, write_ptr);
            ++rr->q;
            rr->new_query = 1;
        }
        else
        {
            trunc_pos = size_to_range(&rr->ix[0], rr->q->beg, space_left);
            if (space_left == buf->alloc
                && less_pair_ordering(&rr->q->beg, &trunc_pos) == 0)
            {
                fprintf(stderr, 
                        "%s: Couldn't fit the line at %Zu:%Zu "
                        "into the allotted space of %Zu\n",
                        __func__, rr->q->beg.hi, rr->q->beg.lo, space_left);
                exit(1);
            }
            
            write_ptr += read_range(&rr->ix[0], rr->q->beg, trunc_pos, write_ptr);
            rr->q->beg = trunc_pos;
            rr->new_query = 0;
        }
        space_left = write_end - write_ptr;
    }
    buf->size = write_ptr - buf->buf;
}
#endif


#define MAX(a, b) ((a) > (b) ? (a) : (b))
#define MIN(a, b) ((a) < (b) ? (a) : (b))

#define MIN_PAIR_ORD(a, b) (less_pair_ordering(&(a), &(b)) < 0 ? (a) : (b))

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
    struct pair_ordering tmp_pos, trunc_pos;

    /* must re-set on each call */
    rr->new_query = 1;
    
    while (rr->new_query && rr->q != rr->qend)
    {
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
                && less_pair_ordering(&rr->q->beg, &trunc_pos) == 0)
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
        space_left = write_end - write_ptr;
    }
    for (i = 0; i != n_ix; ++i)
        bufs[i].size = write_ptr[i] - bufs[i].buf;

    free(span);
    free(write_ptr);
    free(write_end);
}
