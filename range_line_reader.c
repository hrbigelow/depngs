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

    char *write_ptr, *base_end;
    if (*in_alloc == 0)
    {
        *in_alloc = rr->base_size + seen_line_size;
        *in_size = 0;
        *in_buf = malloc(*in_alloc);
        if (! *in_buf)
        {
            fprintf(stderr, "Error: Couldn't allocate %Zu bytes, at %s:%u\n",
                    *in_alloc, __FILE__, __LINE__);
            exit(1);
        }
    }

    write_ptr = *in_buf;
    base_end = write_ptr + rr->base_size;
    
    /* redo the input strategy assuming off_index input */
    while (rr->q != rr->qend)
    {
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

                if (new_line_size != seen_line_size)
                {
                    /* need to realloc */
                    size_t write_pos = write_ptr - *in_buf;
                    seen_line_size = new_line_size;
                    *in_size = rr->base_size + seen_line_size;
                    *in_alloc = *in_size;
                    *in_buf = realloc(*in_buf, *in_alloc);
                    write_ptr = *in_buf + write_pos;
                    base_end = *in_buf + rr->base_size;
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
    }
    if (line_buf)
        free(line_buf);
}
