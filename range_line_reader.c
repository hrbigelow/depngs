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
#define MAX_PAIR_ORD(a, b) (cmp_pair_ordering(&(a), &(b)) < 0 ? (b) : (a))

struct ptr_pair {
    char *cur, *end;
};

struct off_pair {
    off_t beg, end;
};

/* Read as much of several files as possible while reading less than
   max_bytes of any one file.  Start at logical position 'beg' and
   find a new logical position upper bound. Return an array of arrays
   of physical offsets off_ranges[r * n_s + s] = { beg, end }, r is
   range, s is sample. par describes the files and logical ranges from
   those files. pos defines the starting logical position. return the
   end logical position found. */
struct pair_ordering range_line_aux(void *par, unsigned max_bytes,
                                    void *pos,
                                    struct off_pair **off_ranges, 
                                    unsigned *n_off_ranges,
                                    unsigned do_free)
{
    struct range_line_reader_par *rr = par;
    size_t n_ix = rr->n_ix;

    /* advance rr->q and/or rr->q->beg until it coincides with beg */
    struct pair_ordering *pbeg = pos;
    while (rr->qcur != rr->qend 
           && cmp_pair_ordering(&rr->qcur->end, pbeg) == -1) ++rr->qcur;

    off_t max_span;
    off_t *span = malloc(n_ix * sizeof(off_t));

    ptrdiff_t space_left = max_bytes, space_tmp;
    struct pair_ordering tmp_pos, trunc_pos = (struct pair_ordering){ SIZE_MAX, SIZE_MAX };

    unsigned r = 0, i, n_off = n_ix, n_alloc = n_ix;
    struct off_pair *off = malloc(n_off * n_alloc * sizeof(struct off_pair));

    /* must re-set on each call */
    rr->new_query = 1;

    /* initialize this just so that the function returns *pbeg when
       there is no more input to process.  i'm not clear whether this
       makes sense ... */
    struct pair_ordering_range cur_rng = { *pbeg, *pbeg };
    
    while (rr->new_query && rr->qcur != rr->qend)
    {
        cur_rng = (struct pair_ordering_range){ 
            MAX_PAIR_ORD(rr->qcur->beg, *pbeg), rr->qcur->end
        };

        ALLOC_GROW(off, n_off, n_alloc);
        max_span = 0;
        for (i = 0; i != n_ix; ++i)
        {
            span[i] = range_to_size(&rr->ix[i], cur_rng.beg, cur_rng.end);
            max_span = MAX(max_span, span[i]);
        }

        /* fill the buffer as much as possible with the next query range.
           afterwards, qcur now points to the next range to retrieve. */
        unsigned oi;
        if (max_span <= space_left)
        {
            for (i = 0; i != n_ix; ++i)
            {
                oi = r * n_ix + i;
                off[oi].beg = off_lower_bound(&rr->ix[i], cur_rng.beg);
                off[oi].end = off_lower_bound(&rr->ix[i], cur_rng.end);
                space_tmp = off[oi].end - off[oi].beg;
                space_left = MIN(space_left, space_tmp);
            }
            ++rr->qcur;
            rr->new_query = 1;
        }
        else
        {
            for (i = 0; i != n_ix; ++i)
            {
                tmp_pos = size_to_range(&rr->ix[i], cur_rng.beg, space_left);
                trunc_pos = MIN_PAIR_ORD(trunc_pos, tmp_pos);
            }
            if (space_left == max_bytes
                && cmp_pair_ordering(&cur_rng.beg, &trunc_pos) == 0)
            {
                fprintf(stderr, 
                        "%s: Couldn't fit the line at %Zu:%Zu "
                        "into the allotted space of %Zu\n",
                        __func__, cur_rng.beg.hi, cur_rng.beg.lo, space_left);
                exit(1);
            }
            
            for (i = 0; i != n_ix; ++i)
            {
                oi = r * n_ix + i;
                off[oi].beg = off_lower_bound(&rr->ix[i], cur_rng.beg);
                off[oi].end = off_lower_bound(&rr->ix[i], trunc_pos);
                space_tmp = off[oi].end - off[oi].beg;
                space_left = MIN(space_left, space_tmp);
            }
            cur_rng.beg = trunc_pos;
            rr->new_query = 0;
        }
        ++r;
        n_off += n_ix;
    }

    /* free resources.  use max_minus_one to avoid freeing the root */
    struct pair_ordering max_minus_one = { SIZE_MAX, SIZE_MAX - 1 };
    
    if (do_free)
        for (i = 0; i != n_ix; ++i)
            (void)file_bsearch_range_free(&rr->ix[i],
                                          (struct pair_ordering){ 0, 0 },
                                          rr->qcur == rr->qend ? max_minus_one : cur_rng.beg);
    free(span);
    *off_ranges = off;
    *n_off_ranges = r;

    return cur_rng.beg;
}

/* read data starting at beg, using ranges and file handles specified
   in par.  store retrieved data in bufs.  */
void rl_reader(void *par, struct managed_buf *bufs)
{
    struct off_pair *or;
    struct range_line_reader_par *rr = par;
    
    unsigned n_or;
    unsigned max_bytes = bufs[0].alloc;
    (void)range_line_aux(par, max_bytes, &rr->start_pos, &or, &n_or, 1);

    unsigned r, s, ri;
    for (s = 0; s != rr->n_ix; ++s)
    {
        off_t pos = 0;
        for (r = 0; r != n_or; ++r)
        {
            ri = r * rr->n_ix + s;
            off_t bytes = or[ri].end - or[ri].beg;
            fseeko(rr->ix[s].fh, or[ri].beg, SEEK_SET);
            fread(bufs[s].buf + pos, 1, bytes, rr->ix[s].fh);
            pos += bytes;
        }
        bufs[s].size = pos;
    }
    free(or);
}

/* start scanning files at pos, finding the furthest logical position
   end such that the total bytes of any one file in [pos, end) is less
   than max_bytes.  use par to define the total set of valid logical
   ranges, and file handles.  store the original value of pos in par.
   this function should attempt to find these without actually loading
   bulk data, thus saving bandwidth. thread safe.  does not require
   any mutex locking. */
void rl_scanner(void *par, unsigned max_bytes)
{
    struct range_line_reader_par *rr = par;
    struct off_pair *or;
    unsigned n_or;
    rr->end_pos = range_line_aux(par, max_bytes, &rr->start_pos, &or, &n_or, 0);
    free(or);
}


/* transfer thread-local value to global value */
void rl_set_start(void *par, void *pos)
{
    struct range_line_reader_par *rr = par;
    *(struct pair_ordering *)pos = rr->end_pos;
}


/* transfer global value to thread-local value */
void rl_get_start(void *par, void *pos)
{
    struct range_line_reader_par *rr = par;
    rr->start_pos = *(struct pair_ordering *)pos;
}
