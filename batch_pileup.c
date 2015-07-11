/* 
   Populates complete pileup statistics for multiple samples in
   tandem.  Provides access to these statistics for a given sample at
   the current position in three views: basecall counts, (b,q,s)
   triplet counts, or indel counts.  Maintains an internal 'iterator'
   position and provides the user with a 'next' functionality for the
   iterator.

   To use:


*/

#include "compat_util.h"
#include "khash.h"
#include "ksort.h"
#include "cache.h"
#include "bam_reader.h"

#include <stdint.h>
#include <assert.h>

struct indel_seq {
    char is_ins;
    char seq[FLEX_ARRAY]; /* zero-terminated */
};


static inline khint_t indel_hash_func(struct indel_seq *a)
{
    return kh_str_hash_func(a->seq);
}

static inline int indel_hash_equal(struct indel_seq *a, struct indel_seq *b)
{
    return a->is_ins == b->is_ins 
        && kh_str_hash_equal(a->seq, b->seq);
}

struct pbqt {
    uint32_t pos;
    uint32_t tid: 22;
    uint32_t base: 2;
    uint32_t qual: 7;
    uint32_t strand: 1; /* 1 for positive, 0 for negative */
};

union pbqt_key {
    khint64_t k;
    struct pbqt v;
};

struct bqs_count {
    uint32_t base: 2;
    uint32_t qual: 7;
    uint32_t strand: 1;
    uint32_t ct: 22; /* 4,194,304 */
};

struct base_count {
    unsigned ct[4];
};

struct contig_pos {
    unsigned tid, pos;
};

struct pos_base_count {
    struct contig_pos cpos;
    struct base_count bct;
};


struct pos_indel {
    unsigned pos: 32;
    unsigned tid: 22;
    unsigned indel_key: 10; /* 1024 distinct indel types */
};

union pos_key {
    khint64_t k;
    struct contig_pos v;
};

struct indel_count {
    khint_t indel_key;
    unsigned ct;
};

struct indel_count_node {
    struct indel_count_node *next;
    struct indel_count ict;
};

struct pos_indel_count {
    struct contig_pos cpos;
    struct indel_count ict;
};


/* hash type for storing bqt (basecall, quality, strand) counts at
   position.  use  */
KHASH_MAP_INIT_INT64(pbqt_h, unsigned);

/* hash type for storing basecall counts at position */
KHASH_MAP_INIT_INT64(p_h, struct base_count);

/* hash type for storing distinct indel types. */
KHASH_INIT(indel_h, struct indel_seq *, char, 0, indel_hash_func, indel_hash_equal);

/* hash type for tallying indel events. use   */
KHASH_MAP_INIT_INT64(indel_ct_h, struct indel_count_node *);

/* overlapping mates hash, key is q_name.  values are pointers to raw
   BAM records. */
KHASH_MAP_INIT_STR(olap_h, char *);

/* complete tally statistics for the set of BAM records overlapping a
   given set of locus ranges. */
struct tally_stats {
 /* Once tallying is complete, pbqt_hash and indel_ct_hash have
    complete statistics for loci in [tally_beg, tally_end), and
    partial statistics for loci in [tally_end, MAX). */
    khash_t(pbqt_h) *pbqt_hash;
    khash_t(indel_ct_h) *indel_ct_hash;
    /* */
    unsigned n_indel_nodes;

    khash_t(olap_h) *overlap_hash;

    /* summary statistics are compiled for positions in [tally_beg,
       tally_end) */
    khash_t(p_h) *p_hash;
    struct pos_base_count *base_ct, *base_cur, *base_end;
    struct pos_indel_count *indel_ct, *indel_cur, *indel_end;
};


static char **reference_seq;
static unsigned n_reference_seq;

void reference_seq_init(char *fasta_file)
{
}

void reference_seq_free()
{
    unsigned r;
    for (r = 0; r != n_reference_seq; ++r)
        free(reference_seq[r]);
    free(reference_seq);
}

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
            0,
            kh_init(olap_h),
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
        kh_destroy(olap_h, ts->overlap_hash);
        kh_destroy(p_h, ts->p_hash);
        free(ts->base_ct);
        free(ts->indel_ct);
    }
}

static void summarize_base_counts(unsigned s, struct contig_pos tally_end);
static void make_p_array(unsigned s);
static void make_indel_array(unsigned s, struct contig_pos tally_end);
struct contig_pos process_bam_block(char *raw, char *end, struct tally_stats *ts);
static void incr_indel_count(struct tally_stats *ts, struct contig_pos pos,
                             khint_t indel_iter);


static inline int less_contig_pos(struct contig_pos a, struct contig_pos b)
{
    return a.tid < b.tid 
        || (a.tid == b.tid && a.pos < b.pos);
}

static inline int equal_contig_pos(struct contig_pos a, struct contig_pos b)
{
    return a.tid == b.tid && a.pos == b.pos;
}

#define MIN_CONTIG_POS(a, b) (less_contig_pos((a), (b)) ? (a) : (b))
#define MAX_CONTIG_POS(a, b) (less_contig_pos((a), (b)) ? (b) : (a))

/* perform entire tally phase, for basecalls and for indels. bam[s] is
   block of bam records for sample s */
void batch_tally(const struct managed_buf *bam)
{
    unsigned s;
    char *b, *be;

    tls.tally_beg = tls.tally_end;
    struct contig_pos bpos, cpos = { 0, 0 };
    for (s = 0; s != tls.n_samples; ++s)
    {
        b = bam[s].buf;
        be = b + bam[s].size;
        bpos = process_bam_block(b, be, &tls.ts[s]);
        cpos = MAX_CONTIG_POS(cpos, bpos);
    }

    /* if cpos == { 0, 0 } then we are at end of input; summarize all
       remaining statistics */
    if (cpos.tid == 0 && cpos.pos == 0) 
        cpos = (struct contig_pos){ UINT_MAX, UINT_MAX };
        
    tls.tally_end = cpos;
    for (s = 0; s != tls.n_samples; ++s)
    {
        summarize_base_counts(s, tls.tally_end);
        make_p_array(s);
        make_indel_array(s, tls.tally_end);
    }
}


/* provide basecall stats for a sample at current position, or the
   null statistic. */
struct base_count basecall_stats(unsigned s)
{
    static struct base_count null_ct = { { 0, 0, 0, 0 } };
    struct tally_stats *ts = &tls.ts[s];
    if (ts->base_cur == ts->base_end
        || less_contig_pos(tls.base_cur_pos, ts->base_cur->cpos))
        return null_ct;
    else return ts->base_cur->bct;
}


#define MIN_QUAL 0
#define MAX_QUAL 127

/* provide (b,q,s) stats for a sample at current position.  *cts is
   reallocated as necessary. *n_cts set to number of distinct stats
   that are populated. */
void bqs_stats(unsigned s, struct bqs_count **cts, unsigned *n_cts)
{
    unsigned b, q, st;
    union pbqt_key pk;
    khint_t j;
    pk.v.tid = tls.base_cur_pos.tid;
    pk.v.pos = tls.base_cur_pos.pos;
    unsigned n = 0;

    /* should be about 4k.  is this too large? */
    struct bqs_count buf[MAX_QUAL - MIN_QUAL * 4 * 2];

    struct tally_stats *ts = &tls.ts[s];
    for (st = 0; st != 2; ++st)
    {
        pk.v.strand = st;
        for (b = 0; b != 4; ++b)
        {
            pk.v.base = b;
            for (q = MIN_QUAL; q != MAX_QUAL; ++q)
            {
                pk.v.qual = q;
                if ((j = kh_get(pbqt_h, ts->pbqt_hash, pk.k))
                    != kh_end(ts->pbqt_hash))
                    buf[n++] = 
                        (struct bqs_count){ b, q, st, kh_val(ts->pbqt_hash, j) };
            }
        }
    }
    *n_cts = n;
    *cts = realloc(*cts, n * sizeof(struct bqs_count));
    memcpy(*cts, buf, n * sizeof(struct bqs_count));
}


/* provide indel stats for a sample at current position. *cts is
   reallocated as necessary, *n_cts is set to the number of distinct
   indels found. */
void indel_stats(unsigned s, struct indel_count **cts, unsigned *n_cts)
{
    unsigned n = 0;
    struct pos_indel_count *pc = tls.ts[s].indel_cur;

    /* pre-scan to find total number of distinct indels (use linear
       scan since this will be very small) */
    while (pc != tls.ts[s].indel_end &&
           equal_contig_pos(pc->cpos, tls.indel_cur_pos))
        ++n;

    *cts = realloc(*cts, n * sizeof(struct indel_count));
    *n_cts = n;
    unsigned i = 0;
    pc = tls.ts[s].indel_cur;
    while (i != n)
        *cts[i++] = pc++->ict;
}



/* traverse b, tallying the match blocks into ts->pbqt_hash, and the
   indels into ts->indel_ct_hash and tls.indel_hash as necessary */
static void process_bam_stats(bam1_t *b, struct tally_stats *ts)
{
    int32_t tid = b->core.tid, qpos = 0, rpos = b->core.pos, q, r;
    unsigned c;
    uint32_t *cigar = bam_get_cigar(b), op, ln;
    unsigned strand = bam_is_rev(b) ? 0 : 1;
    int ret;

    for (c = 0; c != b->core.n_cigar; ++c)
    {
        op = bam_cigar_op(cigar[c]);
        ln = bam_cigar_oplen(cigar[c]);

        if (op == BAM_CMATCH || op == BAM_CEQUAL || op == BAM_CDIFF)
        {
            int32_t qe = qpos + ln;
            khiter_t it;
            for (q = qpos, r = rpos; q != qe; ++q, ++r)
            {
                union pbqt_key stat = {
                    .v = {
                        r, tid, 
                        bam_seqi(bam_get_seq(b), q), 
                        bam_get_qual(b)[q],
                        strand 
                    }
                };
                if ((it = kh_get(pbqt_h, ts->pbqt_hash, stat.k)) == kh_end(ts->pbqt_hash))
                    it = kh_put(pbqt_h, ts->pbqt_hash, stat.k, &ret);
                kh_val(ts->pbqt_hash, it)++;
            }
        }
        else if (op == BAM_CINS || op == BAM_CDEL)
        {
            /* tally insertion of query. */
            struct indel_seq *isq = malloc(sizeof(struct indel_seq) + ln + 1);
            isq->is_ins = (op == BAM_CINS);
            unsigned i;

            if (isq->is_ins)
                for (i = 0, q = qpos; i != ln; ++i, ++q)
                    isq->seq[i] = seq_nt16_str[bam_seqi(bam_get_seq(b), q)];

            else
                memcpy(isq->seq, reference_seq[tid] + rpos, ln);

            isq->seq[ln] = '\0';

            khiter_t it;
            if ((it = kh_get(indel_h, tls.indel_hash, isq)) == kh_end(tls.indel_hash))
                it = kh_put(indel_h, tls.indel_hash, isq, &ret);
            else free(isq);
            
            struct contig_pos cpos = { tid, rpos };
            incr_indel_count(ts, cpos, it);
        }
        else
            ; /* All other operations (N, S, H, P) do not result in
                 tallying anything */

        /* but, we do advance the query */
        if (bam_cigar_type(op) & 1) qpos += ln;
        if (bam_cigar_type(op) & 2) rpos += ln;
    }
}


/* process an entire BAM block.  maintains a map of possibly
   overlapping reads and resolves their quality scores. returns the
   last position processed, or { 0, 0 } if given empty input (raw ==
   end). */
struct contig_pos process_bam_block(char *raw, char *end, struct tally_stats *ts)
{
    bam1_t b, b_mate;
    int ret;
    char *rtmp;
    struct contig_pos last_pos = { 0, 0 };
    while (raw != end)
    {
        rtmp = bam_parse(raw, &b);
        if (! (b.core.flag & BAM_FPAIRED) || ! (b.core.flag & BAM_FPROPER_PAIR))
        {
            /* either this read is not paired or it is paired but not
               mapped in a proper pair.  tally and don't store. */
            process_bam_stats(&b, ts);
        }
        else
        {
            /* read is paired and mapped in proper pair */
            if (b.core.pos <= b.core.mpos)
            {
                /* b is upstream */
                if (bam_endpos(&b) < b.core.mpos)
                {
                    /* does not overlap. tally and don't store */
                    process_bam_stats(&b, ts);
                }
                else
                {
                    /* b overlaps its mate.  since it is upstream, it
                       shouldn't be in the overlaps hash.  store in
                       overlaps hash. do not tally yet. */
                    khiter_t it = kh_put(olap_h, ts->overlap_hash, bam_get_qname(&b), &ret);
                    kh_val(ts->overlap_hash, it) = rtmp;
                }
            }
            else
            {
                /* b is paired and mapped in proper pair and downstream. the
                   only way we can tell if it overlaps for sure is to search
                   the overlaps hash. */
                khiter_t it =
                    kh_get(olap_h, ts->overlap_hash, bam_get_qname(&b));
                if (it != kh_end(ts->overlap_hash))
                {
                    (void)bam_parse(kh_val(ts->overlap_hash, it), &b_mate);
                    tweak_overlap_quality(&b, &b_mate);
                    process_bam_stats(&b, ts);
                    process_bam_stats(&b_mate, ts);
                    kh_del(olap_h, ts->overlap_hash, it);
                }
                else
                {
                    /* downstream mate does not overlap with its
                       upstream partner */
                    process_bam_stats(&b, ts);
                }
            }
        }
        if (rtmp == end) last_pos = (struct contig_pos){ b.core.tid, b.core.pos };
        raw = rtmp;
    }
    return last_pos;
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
    struct contig_pos kpos;
    for (i = kh_begin(h_src); i != kh_end(h_src); ++i)
    {
        if (! kh_exist(h_src, i)) continue;
        k.k = kh_key(h_src, i);
        kpos = (struct contig_pos){ k.v.tid, k.v.pos };
        if (! less_contig_pos(kpos, tally_end)) continue;

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


static int pos_base_count_less(struct pos_base_count a, struct pos_base_count b)
{
    return less_contig_pos(a.cpos, b.cpos);
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

    khiter_t it;
    unsigned i;
    union pos_key k;
    for (it = kh_begin(ph), i = 0; it != kh_end(ph); ++it)
        if (kh_exist(ph, it))
        {
            k.k = kh_key(ph, it);
            ary[i++] = (struct pos_base_count){ k.v, kh_val(ph, it) };
        }
    
    ks_introsort(pbc_sort, n, ary);
    tls.ts[s].base_cur = ary;
    tls.ts[s].base_end = ary + i;
}

static int pos_indel_count_less(struct pos_indel_count a, struct pos_indel_count b)
{
    return less_contig_pos(a.cpos, b.cpos);
}

KSORT_INIT(pi_sort, struct pos_indel_count, pos_indel_count_less);


/* create a sorted array from sample s's indel_hash */
static void make_indel_array(unsigned s, struct contig_pos tally_end)
{
    khash_t(indel_ct_h) *ih = tls.ts[s].indel_ct_hash;
    unsigned n = tls.ts[s].n_indel_nodes;
    struct pos_indel_count *ary = 
        realloc(tls.ts[s].indel_ct, n * sizeof(struct pos_indel_count));

    khiter_t it;
    unsigned i;
    struct indel_count_node *node;
    union pos_key k;
    for (it = kh_begin(ih), i = 0; it != kh_end(ih); ++it)
        if (kh_exist(ih, it))
        {
            k.k = kh_key(ih, it);
            if (less_contig_pos(k.v, tally_end))
            {
                node = kh_val(ih, it);
                while (node)
                {
                    ary[i++] = (struct pos_indel_count){ k.v, node->ict };
                    node = node->next;
                }
            }
        }

    ks_introsort(pi_sort, n, ary);
    tls.ts[s].indel_ct = ary;
    tls.ts[s].indel_cur = ary;
    tls.ts[s].indel_end = ary + i;
}

/* */
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
static void incr_indel_count(struct tally_stats *ts,
                             struct contig_pos pos,
                             khint_t indel_key)
{
    union pos_key k = { .v = pos };
    khash_t(indel_ct_h) *ih = ts->indel_ct_hash;

    khiter_t i;
    int ret;
    if ((i = kh_get(indel_ct_h, ih, k.k)) == kh_end(ih))
    {
        i = kh_put(indel_ct_h, ih, k.k, &ret);
        kh_val(ih, i) = NULL;
    }

    struct indel_count_node *head, *node;
    head = node = kh_val(ih, i);
    
    /* find node with matching key, or NULL */
    while (node && node->ict.indel_key != indel_key)
        node = node->next;

    if (node) ++node->ict.ct;
    else
    {
        ++ts->n_indel_nodes;
        node = malloc(sizeof(struct indel_count_node));
        *node = (struct indel_count_node){ head, { indel_key, 1 } };
        kh_val(ih, i) = node;
    }
}




/* clear statistics that are no longer needed. call this function
   after next_base() and next_indel() return 1. */
void clear_finished_stats()
{
    khiter_t it;
    unsigned s;
    struct contig_pos kpos;

    for (s = 0; s != tls.n_samples; ++s)
    {
        /* clear finished entries in pbqt hash */
        khash_t(pbqt_h) *ph = tls.ts[s].pbqt_hash;
        union pbqt_key pk;
        for (it = kh_begin(ph); it != kh_end(ph); ++it)
        {
            if (! kh_exist(ph, it)) continue;
            pk.k = kh_key(ph, it);
            kpos = (struct contig_pos){ pk.v.tid, pk.v.pos };
            if (less_contig_pos(kpos, tls.tally_end))
                kh_del(pbqt_h, ph, it);
        }

        /* clear finished entries in indel counts hash */
        khash_t(indel_ct_h) *ih = tls.ts[s].indel_ct_hash;
        union pos_key pk2;
        for (it = kh_begin(ih); it != kh_end(ih); ++it)
        {
            if (! kh_exist(ih, it)) continue;
            pk2.k = kh_key(ih, it);
            if (less_contig_pos(pk2.v, tls.tally_end))
            {
                free_indel_counts(kh_val(ih, it));
                kh_del(indel_ct_h, ih, it);
            }
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
    int next_ ## TYP()                                                  \
    {                                                                   \
        unsigned s;                                                     \
        struct tally_stats *ts;                                         \
                                                                        \
        /* update all TYP_cur pointers */                               \
        struct contig_pos suc =                                         \
            { tls.TYP ## _cur_pos.tid, tls.TYP ## _cur_pos.pos + 1 };   \
        for (s = 0; s != tls.n_samples; ++s)                            \
        {                                                               \
            ts = &tls.ts[s];                                            \
            if (ts->TYP ## _cur != ts->TYP ## _end                      \
                && less_contig_pos(ts->TYP ## _cur->cpos, suc))         \
                ++ts->TYP ## _cur;                                      \
        }                                                               \
        /* recalculate tls.TYP_cur_pos */                               \
        struct contig_pos tmp = { UINT_MAX, UINT_MAX }, tmp2;           \
        for (s = 0; s != tls.n_samples; ++s)                            \
        {                                                               \
            ts = &tls.ts[s];                                            \
            tmp2 = (ts->TYP ## _cur == ts->TYP ## _end)                 \
                ? tmp                                                   \
                : ts->TYP ## _cur->cpos;                                \
            tmp = MIN_CONTIG_POS(tmp, tmp2);                            \
        }                                                               \
        tls.TYP ## _cur_pos = tmp;                                      \
            return                                                      \
                less_contig_pos(tls.TYP ## _cur_pos,                    \
                                (struct contig_pos){ UINT_MAX, UINT_MAX }); \
    }                                                                   \
    

NEXT_AUX(base)

NEXT_AUX(indel)
