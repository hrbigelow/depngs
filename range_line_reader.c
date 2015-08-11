#define _GNU_SOURCE
#define _FILE_OFFSET_BITS 64

#include "range_line_reader.h"
#include "chunk_strategy.h"
#include "virtual_bound.h"

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


extern struct chunk_strategy cs_stats;


struct virt_less_range_par {
    struct pair_ordering_range *ary;
    struct pair_ordering q;
};

/* uses par both as a source of elements elem and a query element q.
   return 1 if q < elem[pos], 0 otherwise.
 */
int less_rng_end(unsigned pos, void *par)
{
    struct virt_less_range_par *vl = par;
    return cmp_pair_ordering(&vl->q, &vl->ary[pos].end) < 0;
}

/* Read as much of several files as possible while reading less than
   max_bytes of any one file.  Start at logical position 'beg' and
   find a new logical position upper bound. Return an array of arrays
   of physical offsets off_ranges[r * n_s + s] = { beg, end }, r is
   range, s is sample. par describes the files and logical ranges from
   those files. pos defines the starting logical position. sets
   cs_stats fields n_bytes_read, n_loci_read. */
void range_line_aux(void *par,
                    unsigned max_bytes,
                    struct off_pair **off_ranges, 
                    unsigned *n_off_ranges,
                    unsigned do_free)
{
    struct range_line_reader_par *rr = par;
    size_t n_ix = rr->n_ix;

    struct pair_ordering_range 
        *qcur, 
        cur_rng = { cs_stats.pos, cs_stats.pos };

    unsigned n_ranges = rr->qend - rr->qbeg;

    /* returns an offset from rr->qbeg such that rr->qbeg[off] is the
       lowest range that extends beyond cs_stats.pos */
    struct virt_less_range_par vpar = { rr->qbeg, cs_stats.pos };
    unsigned q_off =
        virtual_upper_bound(0, n_ranges, less_rng_end, &vpar);
    qcur = rr->qbeg + q_off;

    off_t max_span;
    off_t *span = malloc(n_ix * sizeof(off_t));

    ptrdiff_t space_left = max_bytes;
    struct pair_ordering tmp_pos, trunc_pos = (struct pair_ordering){ SIZE_MAX, SIZE_MAX };

    unsigned r = 0, i, n_off = 0, n_alloc = 0;
    struct off_pair *off = NULL;

    /* must re-set on each call */
    rr->new_query = 1;

    ptrdiff_t max_space_used;

    while (rr->new_query && qcur != rr->qend)
    {
        cur_rng = (struct pair_ordering_range){ 
            MAX_PAIR_ORD(qcur->beg, cs_stats.pos), qcur->end
        };

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
            max_space_used = 0;
            for (i = 0; i != n_ix; ++i)
            {
                oi = r * n_ix + i;
                n_off++;
                ALLOC_GROW(off, n_off, n_alloc);
                off[oi].beg = off_lower_bound(&rr->ix[i], cur_rng.beg);
                off[oi].end = off_lower_bound(&rr->ix[i], cur_rng.end);
                ptrdiff_t used = off[oi].end - off[oi].beg;
                cs_stats.n_bytes_read[i] += used;
                max_space_used = MAX(used, max_space_used);
            }
            assert(max_space_used <= space_left);
            space_left -= max_space_used;
            cs_stats.n_loci_read += qcur->end.lo - cur_rng.beg.lo;
            cur_rng.beg = qcur->end;
            ++qcur;
            rr->new_query = 1;
        }
        else
        {
            /* not enough buffer space to contain the entire current
               query range. determine the furthest possible position
               within the current range that will fit. */
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
            
            max_space_used = 0;
            for (i = 0; i != n_ix; ++i)
            {
                oi = r * n_ix + i;
                n_off++;
                ALLOC_GROW(off, n_off, n_alloc);
                off[oi].beg = off_lower_bound(&rr->ix[i], cur_rng.beg);
                off[oi].end = off_lower_bound(&rr->ix[i], trunc_pos);
                ptrdiff_t used = off[oi].end - off[oi].beg;
                cs_stats.n_bytes_read[i] += used;
                max_space_used = MAX(used, max_space_used);
            }
            assert(max_space_used <= space_left);
            space_left -= max_space_used;
            cs_stats.n_loci_read += trunc_pos.lo - cur_rng.beg.lo;
            cur_rng.beg = trunc_pos;
            rr->new_query = 0;
        }
        ++r;
    }

    /* free resources.  use max_minus_one to avoid freeing the root */
    struct pair_ordering max_minus_one = { SIZE_MAX, SIZE_MAX - 1 };
    
    if (do_free)
        for (i = 0; i != n_ix; ++i)
            (void)file_bsearch_range_free(&rr->ix[i],
                                          (struct pair_ordering){ 0, 0 },
                                          qcur == rr->qend ? max_minus_one : cur_rng.beg);
    free(span);
    *off_ranges = off;
    *n_off_ranges = n_off;

    cs_stats.pos = cur_rng.beg;
}

/* read data starting at beg, using ranges and file handles specified
   in par.  store retrieved data in bufs.  */
void rl_reader(void *par, struct managed_buf *bufs)
{
    struct range_line_reader_par *rr = par;
    
    unsigned r, s, ri;
    for (s = 0; s != rr->n_ix; ++s)
    {
        off_t pos = 0;
        for (r = 0; r != rr->n_offsets / rr->n_ix; ++r)
        {
            ri = r * rr->n_ix + s;
            off_t bytes = rr->offset_pairs[ri].end - rr->offset_pairs[ri].beg;
            fseeko(rr->ix[s].fh, rr->offset_pairs[ri].beg, SEEK_SET);
            fread(bufs[s].buf + pos, 1, bytes, rr->ix[s].fh);
            pos += bytes;
        }
        bufs[s].size = pos;
    }
    free(rr->offset_pairs);
}

/* start scanning files at pos, finding the furthest logical position
   end such that the total bytes of any one file in [pos, end) is less
   than max_bytes.  use par to define the total set of valid logical
   ranges, and file handles.  this function should attempt to find
   these without actually loading bulk data, thus saving
   bandwidth. thread safe.  does not require any mutex locking. */
void rl_scanner(void *par, size_t max_bytes)
{
    struct range_line_reader_par *rr = par;
    size_t bytes_wanted = cs_get_bytes_wanted(rr->n_ix);
    size_t bytes = MIN(bytes_wanted, max_bytes);
    range_line_aux(par, bytes, &rr->offset_pairs, &rr->n_offsets, 0);
}
