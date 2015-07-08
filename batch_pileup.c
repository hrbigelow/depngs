/* 
   Populates complete pileup statistics for multiple samples in
   tandem.  Provides access to these statistics for a given sample at
   the current position in three views: basecall counts, (b,q,s)
   triplet counts, or indel counts.  Maintains an internal 'iterator'
   position and provides the user with a 'next' functionality for the
   iterator.
*/

#include "compat_util.h"
#include "khash.h"
#include "ksort.h"

#include <stdint.h>
#include <assert.h>

struct indel_seq {
    char is_ins;
    char seq[FLEX_ARRAY];
};

struct pbqt {
    unsigned pos: 32;
    unsigned tid: 22;
    unsigned base: 2;
    unsigned qual: 7;
    unsigned strand: 1;
};

union pbqt_key {
    khint64_t k;
    struct pbqt v;
};

struct bqs_count {
    unsigned base: 2;
    unsigned qual: 7;
    unsigned strand: 1;
    unsigned ct;
};

struct base_count {
    unsigned ct[4];
};

struct contig_pos {
    unsigned tid, pos;
};

struct pos_base_count {
    struct contig_pos c;
    struct base_count b;
};


struct pos_indel {
    unsigned pos: 32;
    unsigned tid: 22;
    unsigned indel_key: 10; /* 1024 distinct indel types */
};

union indel_ct_key {
    khint64_t k;
    struct contig_pos v;
};

struct indel_count {
    khint_t indel_key;
    unsigned ct;
};

struct indel_count_node {
    struct indel_count_node *next;
    struct indel_count c;
};

struct pos_indel_count {
    struct contig_pos p;
    struct indel_count c;
};


/* hash type for storing bqt (basecall, quality, strand) counts at
   position */
KHASH_MAP_INIT_INT64(pbqt_h, unsigned);

/* hash type for storing basecall counts at position */
KHASH_MAP_INIT_INT(p_h, struct base_count);

/* hash type for storing distinct indel types */
KHASH_MAP_INIT_INT(indel_h, struct indel_seq *);

/* hash type for tallying indel events. use  */
KHASH_MAP_INIT_INT64(indel_ct_h, struct indel_count_node *);


/* complete tally statistics for the set of BAM records overlapping a
   given set of locus ranges. */
struct tally_stats {
 /* Once tallying is complete, pbqt_hash and indel_ct_hash have
    complete statistics for loci in [tally_beg, tally_end), and
    partial statistics for loci in [tally_end, MAX). */
    khash_t(pbqt_h) *pbqt_hash;
    khash_t(indel_ct_h) *indel_ct_hash;

    /* summary statistics are compiled for positions in [tally_beg,
       tally_end) */
    khash_t(p_h) *p_hash;
    struct pos_base_count *base_ct, *base_cur, *base_end;
    struct pos_indel_count *indel_ct, *indel_cur, *indel_end;
};

static __thread struct tls {
    khash_t(indel_h) *indel_hash;
    unsigned n_samples;

    struct contig_pos tally_beg, tally_end;
    struct tally_stats *ts; /* tally stats for each sample */
    struct contig_pos base_cur_pos; /* current base position being
                                       queried. Initialized in ? */
    struct contig_pos indel_cur_pos; /* current indel position being queried */
} tls;


void batch_pileup_init(unsigned n_samples)
{
    tls.n_samples = n_samples;
    tls.ts = malloc(n_samples * sizeof(struct tally_stats));
    unsigned s;
    for (s = 0; s != n_samples; ++s)
        tls.ts[s] = (struct tally_stats){
            kh_init(pbqt_h),
            kh_init(indel_ct_h),
            kh_init(p_h),
            NULL, NULL, NULL,
            NULL, NULL, NULL
        };
    tls.tally_end = (struct contig_pos){ 0, 0 };
}


void batch_pileup_free()
{
    unsigned s;
    for (s = 0; s != tls.n_samples; ++s)
    {
        struct tally_stats *ts = &tls.ts[s];
        kh_destroy(pbqt_h, ts->pbqt_hash);
        kh_destroy(indel_ct_h, ts->indel_ct_hash);
        kh_destroy(p_h, ts->p_hash);
        free(ts->base_ct);
        free(ts->indel_ct);
    }
}

static void summarize_base_counts(unsigned s, struct contig_pos tally_end);
static void make_p_array(unsigned s, struct contig_pos tally_end);
static void make_indel_array(unsigned s, struct contig_pos tally_end);


/* perform entire tally phase. bam[s] is block of bam records for
   sample s */
void tally(const struct managed_buf *bam)
{
    unsigned s;
    char *b, *b_pre, *b_end;
    bam_rec_t rec;

    tls.tally_beg = tls.tally_end;
    struct contig_pos cpos = { 0, 0 };
    for (s = 0; s != tls.n_samples; ++s)
    {
        b = bam[s].buf;
        b_end = bam[s].buf + bam[s].size;
        while (b != b_end)
        {
            b_pre = b;
            b = tally_bam(s, b);
        }
        cpos = MAX_POS(cpos, bam_pos(b_pre));
    }
    /* if pos == { 0, 0 } then we are at end of input; summarize all
       remaining statistics */
    if (cpos.tid == 0 && cpos.pos == 0) 
        cpos = (struct contig_pos){ UINT_MAX, UINT_MAX };
        
    tls.tally_end = pos;
    for (s = 0; s != tls.n_samples; ++s)
    {
        summarize_base_counts(s, tls.tally_end);
        make_p_array(s, tls.tally_end);
        make_indel_array(s, tls.tally_end);
    }
}


/* provide basecall stats for a sample at current position, or the
   null statistic. */
struct base_count basecall_stats(unsigned s)
{
    static struct base_count null_ct = { 0, 0, 0, 0 };
    struct tally_stats *ts = &tls.ts[s];
    if (ts->base_cur == ts->base_end
        || less_contig_pos(tls.base_cur_pos, ts->base_cur->pos))
        return null_ct;
    else return ts->base_cur->b;
}


#define MIN_QUAL 0
#define MAX_QUAL 127

/* provide (b,q,s) stats for a sample at current position.  cts must
   have at least 2 * 4 * N_QUAL */
void bqs_stats(unsigned s, struct bqs_count *cts, unsigned *n_cts)
{
    unsigned b, q, s;
    union pbqt_key pk;
    khint_t j;
    pk.v.tid = tls.base_cur_pos.tid;
    pk.v.pos = tls.base_cur_pos.pos;
    unsigned n = 0;
    struct tally_stats *ts = &tls.ts[s];
    for (s = 0; s != 2; ++s)
    {
        pk.v.strand = s;
        for (b = 0; b != 4; ++b)
        {
            pk.v.base = b;
            for (q = MIN_QUAL; q != MAX_QUAL; ++q)
            {
                pk.v.qual = q;
                if ((j = kh_get(pbqt_h, ts->pbqt_hash, pk.k))
                    != kh_end(ts->pbqt_hash))
                    cts[n++] = 
                        (struct bqs_count){ b, q, s, kh_val(ts->pbqt_hash, j) };
            }
        }
    }
    *n_cts = n;
}


/* provide indel stats for a sample at current position */
void indel_stats(unsigned s)
{
    
}

/* tally one BAM record for sample s, updating tls hashes. return
   pointer to next BAM record. */
static char *tally_bam(unsigned s, char *bam)
{
    /* traverse the BAM record */
    
}


/* condense to p_hash */
static void summarize_base_counts(unsigned s, struct contig_pos tally_end)
{
    union pbqt_key k;
    unsigned ct;
    khash_t(pbqt_h) *h_src = tls.ts[s].pbqt_hash;
    khash_t(p_h) *h_trg = tls.ts[s].p_hash;
    kh_clear(p_h, h_trg);
    khint_t i, j;
    int ret;
    for (i = kh_begin(h_src); i != kh_end(h_src); ++i)
    {
        if (! kh_exist(h_src, i)) continue;
        k.k = kh_key(h_src, i);
        if (k.v.pos >= tally_end) continue;

        ct = kh_val(h_src, i);
        j = kh_get(p_h, h_trg, k.v.pos);
        if (j == kh_end(h_trg)) 
        {
            j = kh_put(p_h, h_trg, k.v.pos, &ret);
            assert(ret == 0);
            kh_val(h_trg, j) = (struct base_count){ { 0, 0, 0, 0 } };
        }
        kh_val(h_trg, j).ct[k.v.base] += ct;
    }
}


static inline int contig_pos_less(struct contig_pos a, struct contig_pos b)
{
    return a.tid < b.tid 
        || (a.tid == b.tid && a.pos < b.pos);
}

static int pos_base_count_less(struct pos_base_count a, struct pos_base_count b)
{
    return contig_pos_less(a.c, b.c);
}

KSORT_INIT(pbc_sort, struct pos_base_count, pos_base_count_less);

/* create a sorted array from sample s's p_hash.  Note: p_hash is only
   summarized for [tally_beg, tally_end). */
static void make_p_array(unsigned s)
{
    khash_t(p_h) *ph = tls.ts[s].p_hash;
    unsigned n = kh_size(ph);
    struct pos_base_count *ary = 
        realloc(tls.ts[s].base_ct, n * sizeof(struct pos_base_count));

    khint_t k;
    unsigned i;
    for (k = kh_begin(ph), i = 0; k != kh_end(ph); ++k)
        if (kh_exist(ph, k))
            ary[i++] = (struct pos_base_count){ kh_key(ph, k), kh_val(ph, k) };
    
    ks_introsort(pbc_sort, n, ary);
    tls.ts[s].base_cur = ary;
    tls.ts[s].base_end = ary + i;
}

static int pos_indel_count_less(struct pos_indel_count a, struct pos_indel_count b)
{
    return contig_pos_less(a.c, b.c);
}

KSORT_INIT(pi_sort, struct pos_indel_count, pos_indel_count_less);

/* create a sorted array from sample s's indel_hash */
static void make_indel_array(unsigned s, struct contig_pos tally_end)
{
    khash_t(indel_h) *ih = tls.indel_hash;
    unsigned n = kh_size(ih);
    struct pos_indel_count *ary = 
        realloc(tls.ts[s].indel_ct, n * sizeof(struct pos_indel_count));

    khint_t k;
    unsigned i;
    for (k = kh_begin(ih), i = 0; k != kh_end(ih); ++k)
        if (kh_exist(ih, k) && kh_key(ih, k).v.pos < tally_end)
            ary[i++] = (struct pos_indel_count){ kh_key(ih, k), kh_val(ih, k) };

    ks_introsort(pi_sort, n, ary);
    tls.ts[s].indel_ct = ary;
    tls.ts[s].indel_cur = ary;
    tls.ts[s].indel_end = ary + i;
}


static void free_indel_counts(struct indel_count_node *nd)
{
    assert(nd != NULL);
    struct indel_count_node *t = nd->next;
    do {
        free(nd);
        nd = t;
        t = t->next;
    } while (t);
}


/* increment a count of a particular indel, represented by indel_key,
   in a sample s at pos. */
static void incr_indel_count(unsigned s, struct contig_pos pos, khint_t indel_key)
{
    struct tally_stats *ts = &tls.ts[s];
    khiter_t i;
    union indel_ct_key k = { .v = pos };
    if ((i = kh_get(indel_ct_h, ts->indel_ct_hash, k.k)) == kh_end(ts->indel_ct_hash))
        i = kh_put(indel_ct_h, ts->indel_ct_hash, k.k);
    kh_val(ts->indel_ct_hash, i) =
}

/* clear statistics that are no longer needed */
static void clear_finished_stats(struct contig_pos tally_end)
{
    khiter_t i;
    khint_t pos;
    unsigned s;
    union pbqt_key pk;
    
    for (s = 0; s != tls.n_samples; ++s)
    {
        /* clear finished entries in pbqt hash */
        khash_t(pbqt_h) *h = tls.ts[s].pbqt_hash;
        for (i = kh_begin(h); i != kh_end(h); ++i)
        {
            if (! kh_exist(h, i)) continue;
            pk.k = kh_key(h, i);
            if (pk.v.pos < tally_end) kh_del(h, i);
        }

        /* clear finished entries in indel counts hash */
        khash_t(indel_ct_h) *ih = tls.ts[s].indel_ct_hash;
        for (i = kh_begin(ih); i != kh_end(ih); ++i)
        {
            if (! kh_exist(ih, i)) continue;
            pos = kh_key(ih, i);
            if (pos < tally_end) kh_del(ih, i);
        }

        /* completely clear temporary data */
        kh_clear(p_h, tls.ts[s].p_hash);
        free(tls.ts[s].base_ct);
        free(tls.ts[s].indel_ct);
    }
}

/* advance TYP (base or indel) iterator position. return 1 if reached
   the end. */
#define NEXT_AUX(TYP)                                                   \
    int next_##TYP##()                                                  \
    {                                                                   \
        unsigned s;                                                     \
        struct tally_stats *ts;                                         \
                                                                        \
        /* update all TYP_cur pointers */                               \
        struct contig_pos suc =                                         \
            { tls.##TYP##_cur_pos.tid, tls.##TYP##_cur_pos.pos + 1 };   \
        for (s = 0; s != tls.n_samples; ++s)                            \
        {                                                               \
            ts = &tls.ts[s];                                            \
            if (ts->##TYP##_cur != ts->##TYP##_end                      \
                && less_contig_pos(ts->##TYP##_cur, suc))               \
                ++ts->##TYP##_cur;                                      \
        }                                                               \
        /* recalculate tls.##TYP##_cur_pos */                           \
        struct contig_pos tmp = { UINT_MAX, UINT_MAX };                 \
        for (s = 0; s != tls.n_samples; ++s)                            \
        {                                                               \
            ts = &tls.ts[s];                                            \
            tmp = min_contig_pos(tmp, ts->##TYP##_cur == ts->##TYP##_end \
                                 ? tmp : *ts->##TYP##_cur);             \
        }                                                               \
        tls.##TYP##_cur_pos = tmp;                                      \
            return less_contig_pos(tls.##TYP##_cur_pos,                 \
                                   (struct contig_pos){ UINT_MAX, UINT_MAX }); \
    }                                                                   \
    

NEXT_AUX(base)

NEXT_AUX(indel)
