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

/* tallys statistics that allow the reader to adjust its strategy
   during program execution. one of two strategies is used:

   1.  if do_range_estimation is set, we have been given a list of
       locus ranges to process.  n_loci_total is initialized at the
       start. n_loci_left, n_bytes_read are maintained.  number of
       bytes left is estimated from these in rl_scanner.

   2.  otherwise, we are processing the whole file.  n_bytes_total is
       initialized at the start.  n_bytes_left is maintained directly.
       n_loci_total and n_loci_left are not used.

 */

static struct {
    /* marker that informs all threads where to resume reading */
    struct pair_ordering pos;

    unsigned do_range_estimation;

    unsigned long *n_bytes_total;
    unsigned long *n_bytes_read;
    unsigned long n_loci_total;
    unsigned long n_loci_left; /* updated by range_line_aux */
} reader_state;


/* Read as much of several files as possible while reading less than
   max_bytes of any one file.  Start at logical position 'beg' and
   find a new logical position upper bound. Return an array of arrays
   of physical offsets off_ranges[r * n_s + s] = { beg, end }, r is
   range, s is sample. par describes the files and logical ranges from
   those files. pos defines the starting logical position. sets
   reader_state fields n_bytes_read, n_loci_left. */
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
        cur_rng = { reader_state.pos, reader_state.pos };

    unsigned n_ranges = rr->qend - rr->qbeg;

    /* find the range that contains cur_rng */
    qcur = bsearch(&cur_rng, rr->qbeg, n_ranges, sizeof(*rr->qbeg), 
                   cmp_pair_ordering_contained);
    if (! qcur) qcur = rr->qend;

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
            MAX_PAIR_ORD(qcur->beg, reader_state.pos), qcur->end
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
                reader_state.n_bytes_read[i] += used;
                max_space_used = MAX(used, max_space_used);
            }
            assert(max_space_used <= space_left);
            space_left -= max_space_used;
            reader_state.n_loci_left -= qcur->end.lo - cur_rng.beg.lo;
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
                reader_state.n_bytes_read[i] += used;
                max_space_used = MAX(used, max_space_used);
            }
            assert(max_space_used <= space_left);
            space_left -= max_space_used;
            reader_state.n_loci_left -= trunc_pos.lo - cur_rng.beg.lo;
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

    reader_state.pos = cur_rng.beg;
}

/* call this if a range file is given */
void rl_init_by_range(unsigned n_loci_total, unsigned n_files)
{
    reader_state.pos = (struct pair_ordering){ 0, 0 };
    reader_state.do_range_estimation = 1;
    reader_state.n_bytes_total = NULL;
    reader_state.n_bytes_read = calloc(n_files, sizeof(reader_state.n_bytes_read[0]));
    reader_state.n_loci_total = n_loci_total;
    reader_state.n_loci_left = n_loci_total;
}

/* call this if no range file is given */
void rl_init_whole_file(struct file_bsearch_index *ix, unsigned n_files)
{
    reader_state.pos = (struct pair_ordering){ 0, 0 };
    reader_state.do_range_estimation = 0;
    reader_state.n_bytes_total = malloc(n_files * sizeof(reader_state.n_bytes_total[0]));
    reader_state.n_bytes_read = calloc(n_files, sizeof(reader_state.n_bytes_read[0]));
    reader_state.n_loci_total = 0;
    reader_state.n_loci_left = 0;

    unsigned f;
    for (f = 0; f != n_files; ++f)
        reader_state.n_bytes_total[f] = 
            ix[f].root->end_offset - ix[f].root->start_offset;
}

void rl_free()
{
    if (reader_state.n_bytes_total) free(reader_state.n_bytes_total);
    free(reader_state.n_bytes_read);
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
   ranges, and file handles.  store the original value of pos in par.
   this function should attempt to find these without actually loading
   bulk data, thus saving bandwidth. thread safe.  does not require
   any mutex locking. */

/* This size is big enough for a pileup line and small enough to be a
   quickly-processed chunk, increasing the chance that multiple
   threads can keep working until the end of the input. */

/* See comments on struct reader_state for description of chunking
   strategy. */

/* If we have less than 1GB of input to go, switch to using small
   chunks of SMALL_CHUNK size */
#define MAX_BYTES_SMALL_CHUNK 1e9
#define SMALL_CHUNK 5e6

/* If in do_range_estimation mode, and we don't yet have an estimate
   for n_bytes_per_locus, use this as a default value. */
#define DEFAULT_BYTES_PER_LOCUS 100
void rl_scanner(void *par, unsigned max_bytes)
{
    struct range_line_reader_par *rr = par;
    unsigned f;
    unsigned long most_bytes_left = 0;
    if (reader_state.do_range_estimation)
    {
        unsigned long loci_read = reader_state.n_loci_total - reader_state.n_loci_left;
        for (f = 0; f != rr->n_ix; ++f)
        {
            unsigned n_bytes_per_locus = 
                loci_read == 0
                ? DEFAULT_BYTES_PER_LOCUS
                : reader_state.n_bytes_read[f] / loci_read;
            
            most_bytes_left = MAX(reader_state.n_loci_left * n_bytes_per_locus, 
                                  most_bytes_left);
        }
    }
    else
        for (f = 0; f != rr->n_ix; ++f)
            most_bytes_left = MAX(most_bytes_left,
                                  reader_state.n_bytes_total[f]
                                  - reader_state.n_bytes_read[f]);

    unsigned bytes_wanted = most_bytes_left > MAX_BYTES_SMALL_CHUNK ? max_bytes : SMALL_CHUNK;
    range_line_aux(par, bytes_wanted, &rr->offset_pairs, &rr->n_offsets, 0);
}
