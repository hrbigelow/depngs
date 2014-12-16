#define _GNU_SOURCE
#define _FILE_OFFSET_BITS 64

#include "range_line_reader.h"

#include <stddef.h>
#include <stdlib.h>
#include <assert.h>
#include <string.h>
#include <stdio.h>

/* Conforms to the thread_queue_reader_t expectation. Manages
   in_buf/in_size/in_alloc.  Assumes line-based input. */
void range_line_reader(void *par, char **in_buf, size_t *in_size, size_t *in_alloc)
{
    struct range_line_reader_par *rr = par;
    size_t new_line_size, seen_line_size = rr->max_line_size;
    off_t beg_off, end_off, span;
    ptrdiff_t base_left;
    char *line_buf = NULL;

    assert(*in_alloc > 0);

    size_t end_pos = *in_alloc - 
        (*in_alloc > rr->max_line_size ? rr->max_line_size : 0);

    char *write_ptr = *in_buf, *base_end = *in_buf + end_pos;

    rr->new_query = 1; /* on each call, we must re-set this, since
                          beg_off and end_off values are not
                          maintained between calls */

    /* redo the input strategy assuming off_index input */
    while ((base_left = base_end - write_ptr) > 0 && rr->q != rr->qend)
    {
        /* find file offsets for current query, creating index
           nodes and updating ix in the process */
        if (rr->new_query)
        {
            rr->ix = find_loose_index(rr->ix, rr->q->beg, rr->fh);
            beg_off = off_lower_bound(rr->ix, rr->q->beg);
            rr->ix = find_loose_index(rr->ix, rr->q->end, rr->fh);
            end_off = off_upper_bound(rr->ix, rr->q->end);

            fseeko(rr->fh, beg_off, SEEK_SET);
        }

        /* fill the buffer as much as possible with the next query range.
           afterwards, q now points to the next range to retrieve. */
        if ((span = end_off - beg_off) < base_left)
        {
            write_ptr += fread(write_ptr, 1, span, rr->fh);
            ++rr->q;
            rr->new_query = 1;
        }
        else
        {
            /* partially consume the query range.  read up to the
               base buffer, then read the next line fragment.
               realloc both in_buf and line_buf as necessary */
            write_ptr += fread(write_ptr, 1, base_left, rr->fh);

            ssize_t line_length = getline(&line_buf, &new_line_size, rr->fh);
            assert(line_length >= 0);

            if (new_line_size > seen_line_size)
            {
                /* need to realloc */
                size_t write_pos = write_ptr - *in_buf;
                *in_alloc += (new_line_size - seen_line_size);
                *in_buf = realloc(*in_buf, *in_alloc);
                write_ptr = *in_buf + write_pos;
                base_end = *in_buf + *in_alloc - new_line_size;
                seen_line_size = new_line_size;
            }
            strcpy(write_ptr, line_buf);
            write_ptr += line_length;
            beg_off = ftello(rr->fh);

            /* update q->beg to the position just after the last line read */
            char *last_line = memrchr(*in_buf, '\n', write_ptr - *in_buf - 1) + 1;
            rr->q->beg = rr->get_line_ord(last_line);
            ++rr->q->beg.lo; /* lo is the locus position, hi is the contig */
            rr->new_query = 0;
        }
    }
    if (line_buf)
        free(line_buf);
    *in_size = write_ptr - *in_buf;
}
