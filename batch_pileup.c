/* 
   Populates complete pileup statistics for multiple samples in
   tandem.  Provides access to these statistics for a given sample at
   the current position in three views: basecall counts, (b,q,s)
   triplet counts, or indel counts.  Maintains an internal 'iterator'
   position and provides the user with a 'next' functionality for the
   iterator.
*/

#include "batch_pileup.h"
#include "htslib/faidx.h"
#include "ksort.h"
#include "cache.h"
#include "bam_reader.h"

#include <stdint.h>
#include <assert.h>
#include <ctype.h>


static inline khint_t
indel_hash_func(struct indel_seq *a)
{
    return kh_str_hash_func(a->seq);
}

static inline int
indel_hash_equal(struct indel_seq *a, struct indel_seq *b)
{
    return a->is_ins == b->is_ins 
        && kh_str_hash_equal(a->seq, b->seq);
}

struct pbqt {
    uint32_t pos;
    uint32_t tid: 20; /* 1,048,576 */
    uint32_t base: 4; /* BAM encoding. use hts.h: seq_nt16_str[base]
                         to get the base */
    uint32_t qual: 7; /* numeric quality score */
    uint32_t strand: 1; /* 1 for positive, 0 for negative */
};

union pbqt_key {
    khint64_t k;
    struct pbqt v;
};

struct contig_pos {
    unsigned tid, pos;
};

struct pos_base_count {
    struct contig_pos cpos;
    struct base_count bct;
};


struct pos_bqs_count {
    struct contig_pos cpos;
    struct bqs_count bqs_ct;
};


struct pos_indel {
    unsigned pos: 32;
    unsigned tid: 22;
    unsigned indel_itr: 10; /* 1024 distinct indel types */
};

union pos_key {
    khint64_t k;
    struct contig_pos v;
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
   position.  uses 'union pbqt_key' as the key  */
KHASH_MAP_INIT_INT64(pbqt_h, unsigned);

/* hash type for storing basecall counts at position. uses 'union
   pos_key' as the key. */
KHASH_MAP_INIT_INT64(pb_h, struct base_count);

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
       complete statistics for loci < tally_end, and
       partial statistics for loci in [tally_end, MAX). */
    khash_t(pbqt_h) *pbqt_hash;
    khash_t(indel_ct_h) *indel_ct_hash;
    /* */
    unsigned n_indel_nodes;
    khash_t(olap_h) *overlap_hash;

    /* summary statistics are compiled for positions < tally_end */
    struct pos_base_count *base_ct, *base_cur, *base_end;
    struct pos_bqs_count *bqs_ct, *bqs_cur, *bqs_end;
    struct pos_indel_count *indel_ct, *indel_cur, *indel_end;
};



static __thread struct tls {
    khash_t(indel_h) *indel_hash;
    unsigned n_samples;
    struct contig_pos tally_end;
    struct tally_stats *ts; /* tally stats for each sample */
    struct contig_pos cur_pos; /* current position being queried */
} tls;



static unsigned min_quality_score;


void
batch_pileup_init(unsigned min_qual, const char *fasta_file)
{
    min_quality_score = min_qual;
    unsigned do_load_seqs = 1;
    genome_init(fasta_file, do_load_seqs);
}


void
batch_pileup_free()
{
    genome_free();
}


#define TID_UNSET UINT_MAX

/* */
void
batch_pileup_thread_init(unsigned n_samples)
{
    tls.indel_hash = kh_init(indel_h);
    tls.n_samples = n_samples;
    tls.tally_end = (struct contig_pos){ 0, 0 };
    tls.cur_pos = (struct contig_pos){ TID_UNSET, 0 };
    tls.ts = malloc(n_samples * sizeof(struct tally_stats));
    unsigned s;
    for (s = 0; s != n_samples; ++s)
        tls.ts[s] = (struct tally_stats){
            .pbqt_hash = kh_init(pbqt_h),
            .indel_ct_hash = kh_init(indel_ct_h),
            .n_indel_nodes = 0,
            .overlap_hash = kh_init(olap_h),
            .base_ct = NULL, 
            .base_cur = NULL, 
            .base_end = NULL,
            .bqs_ct = NULL,
            .bqs_cur = NULL,
            .bqs_end = NULL,
            .indel_ct = NULL, 
            .indel_cur = NULL, 
            .indel_end = NULL
            
        };
}


void
batch_pileup_thread_free()
{
    khiter_t it;
    for (it = kh_begin(tls.indel_hash); it != kh_end(tls.indel_hash); ++it)
        if (kh_exist(tls.indel_hash, it))
            free(kh_key(tls.indel_hash, it));

    kh_destroy(indel_h, tls.indel_hash);

    unsigned s;
    for (s = 0; s != tls.n_samples; ++s) {
        struct tally_stats *ts = &tls.ts[s];
        kh_destroy(pbqt_h, ts->pbqt_hash);
        kh_destroy(indel_ct_h, ts->indel_ct_hash);
        kh_destroy(olap_h, ts->overlap_hash);
        free(ts->base_ct);
        free(ts->bqs_ct);
        free(ts->indel_ct);
    }
}

struct contig_pos
process_bam_block(char *raw, char *end, struct tally_stats *ts);

static void
incr_indel_count(struct tally_stats *ts, struct contig_pos pos,
                 khiter_t indel_itr);


static inline int
less_contig_pos(struct contig_pos a, struct contig_pos b)
{
    return a.tid < b.tid 
        || (a.tid == b.tid && a.pos < b.pos);
}

static inline int
equal_contig_pos(struct contig_pos a, struct contig_pos b)
{
    return a.tid == b.tid && a.pos == b.pos;
}

#define MIN_CONTIG_POS(a, b) (less_contig_pos((a), (b)) ? (a) : (b))
#define MAX_CONTIG_POS(a, b) (less_contig_pos((a), (b)) ? (b) : (a))

/* perform entire tally phase, for basecalls and for indels for one
   sample. update tls.tally_end to reflect furthest start position of
   any bam record seen. */
void
pileup_tally_stats(const struct managed_buf bam, unsigned s)
{
    char *b = bam.buf, *be = b + bam.size;
    struct contig_pos bpos = 
        process_bam_block(b, be, &tls.ts[s]);
    tls.tally_end = MAX_CONTIG_POS(tls.tally_end, bpos);
}



/* return BAM bitflag-encoded form of current refbase. */
static unsigned
get_cur_refbase_code16()
{
    char refbase = reference_seq.contig[tls.cur_pos.tid].seq[tls.cur_pos.pos];
    return seq_nt16_table[(int)refbase];
}


/* */
static unsigned
get_cur_refbase_code5()
{
    return seq_nt16_int[get_cur_refbase_code16()];
}

/* provide basecall stats for a sample at current position, or the
   null statistic. */
struct base_count
pileup_current_basecalls(unsigned s)
{
    static struct base_count null_ct = { { 0, 0, 0, 0 }, 0, 0 };
    static struct base_count refsam_ct[] = {
        { { 1e6, 0, 0, 0 }, 0, 1e6 },
        { { 0, 1e6, 0, 0 }, 0, 1e6 },
        { { 0, 0, 1e6, 0 }, 0, 1e6 },
        { { 0, 0, 0, 1e6 }, 0, 1e6 },
        { { 0, 0, 0, 0 }, 0, 0 }
    };

    if (s == REFERENCE_SAMPLE)
        return refsam_ct[get_cur_refbase_code5()];

    struct tally_stats *ts = &tls.ts[s];
    if (ts->base_cur == ts->base_end
        || less_contig_pos(tls.cur_pos, ts->base_cur->cpos))
        return null_ct;

    else
        return ts->base_cur->bct;
}


static int
pos_bqs_count_less(struct pos_bqs_count a, struct pos_bqs_count b)
{
    return less_contig_pos(a.cpos, b.cpos);
}

KSORT_INIT(pbqt_sort, struct pos_bqs_count, pos_bqs_count_less);


/* transform contents of pbqt hash into a position-sorted array of bqt
   counts. */
void
pileup_prepare_bqs(unsigned s)
{
    struct tally_stats *ts = &tls.ts[s];
    khash_t(pbqt_h) *ph = ts->pbqt_hash;
    unsigned n = kh_size(ph);
    ts->bqs_ct = realloc(ts->bqs_ct, n * sizeof(struct pos_bqs_count));

    khiter_t it;
    unsigned i;
    union pbqt_key pk;
    unsigned ct;
    struct contig_pos pos;
    for (it = kh_begin(ph), i = 0; it != kh_end(ph); ++it) {
        if (! kh_exist(ph, it)) continue;
        pos = (struct contig_pos){ pk.v.tid, pk.v.pos };
        if (less_contig_pos(pos, tls.tally_end)) {
            pk.k = kh_key(ph, it);
            ct = kh_val(ph, it);
            ts->bqs_ct[i++] = (struct pos_bqs_count){
                .cpos = { pk.v.tid, pk.v.pos },
                .bqs_ct = { .base = pk.v.base,
                            .qual = pk.v.qual,
                            .strand = pk.v.strand,
                            .ct = ct }
            };
        }
    }
    ks_introsort(pbqt_sort, i, ts->bqs_ct);
    ts->bqs_cur = ts->bqs_ct;
    ts->bqs_end = ts->bqs_ct + i;
}


/* provide (b,q,s) stats for a sample at current position.  *cts is
   reallocated as necessary. *n_cts set to number of distinct stats
   that are populated. complete statistics are provided; there is no
   filtering by quality score. */
void
pileup_current_bqs(unsigned s, struct bqs_count **cts, unsigned *n_cts)
{
    static struct bqs_count ref_bqs[] = {
        { 0, 50, 0, 1e6 },
        { 1, 50, 0, 1e6 },
        { 2, 50, 0, 1e6 },
        { 3, 50, 0, 1e6 },
        { 15, 0, 0, 0 } /* zero counts of base 'N' */
    };

    if (s == REFERENCE_SAMPLE) {
        *n_cts = 1;
        *cts = realloc(*cts, 1 * sizeof(struct bqs_count));
        (*cts)[0] = ref_bqs[get_cur_refbase_code5()];
        return;
    }

    struct tally_stats *ts = &tls.ts[s];
    struct pos_bqs_count *pc = ts->bqs_cur;
    while (pc != ts->bqs_end && 
           equal_contig_pos(pc->cpos, tls.cur_pos))
        ++pc;
    unsigned n = pc - ts->bqs_cur;
        
    *n_cts = n;
    *cts = realloc(*cts, n * sizeof(struct bqs_count));
    pc = ts->bqs_cur;
    unsigned i;
    for (i = 0; i != n; ++i, ++pc)
        (*cts)[i] = pc->bqs_ct;
}


/* provide indel stats for a sample at current position. *cts is
   reallocated as necessary, *n_cts is set to the number of distinct
   indels found. */
void
pileup_current_indels(unsigned s, struct indel_count **cts, unsigned *n_cts)
{
    if (s == REFERENCE_SAMPLE) {
        *n_cts = 0;
        return;
    }
    
    /* pre-scan to find total number of distinct indels (use linear
       scan since this will be very small) */
    struct pos_indel_count *pc = tls.ts[s].indel_cur;
    while (pc != tls.ts[s].indel_end &&
           equal_contig_pos(pc->cpos, tls.cur_pos))
        ++pc;
    unsigned n = pc - tls.ts[s].indel_cur;

    *cts = realloc(*cts, n * sizeof(struct indel_count));
    *n_cts = n;
    unsigned i = 0;
    pc = tls.ts[s].indel_cur;
    while (i != n)
        (*cts)[i++] = pc++->ict;
}


/* merge information in indel counts 1 and 2, producing counts of
   pairs based on indel type. */
void
pileup_make_indel_pairs(struct indel_count *cts1, unsigned n_cts1,
                        struct indel_count *cts2, unsigned n_cts2,
                        struct indel_pair_count **pair_cts, unsigned *n_pair_cts)
{
    unsigned max_n_events = n_cts1 + n_cts2 + 2;
    *pair_cts = realloc(*pair_cts, max_n_events * sizeof(struct indel_pair_count));

    struct indel_count
        *ic0 = cts1, *ie0 = ic0 + n_cts1,
        *ic1 = cts2, *ie1 = ic1 + n_cts2;

    struct indel_pair_count *ip = *pair_cts;

    /* */
    while (ic0 != ie0 || ic1 != ie1) {
        ip->count[0] = ic0 != ie0 
            && (ic1 == ie1 || ic0->indel_itr <= ic1->indel_itr) ? ic0->ct : 0;
        
        ip->count[1] = ic1 != ie1
            && (ic0 == ie0 || ic1->indel_itr <= ic0->indel_itr) ? ic1->ct : 0;
        
        if (ip->count[0] != 0) { ip->indel_id = ic0->indel_itr; ++ic0; }
        if (ip->count[1] != 0) { ip->indel_id = ic1->indel_itr; ++ic1; }
        ++ip;
    }            
    *n_pair_cts = ip - *pair_cts;
}


/* */
static inline char
bqs_count_to_call(struct bqs_count bc, unsigned refbase_i)
{
    static char match[] = ".,";
    return bc.base == refbase_i
        ? match[bc.strand]
        : seq_nt16_str[bc.base];
}

/* convert a BAM quality to a character */
static inline char
qual_to_char(unsigned q)
{
    return (char)(q + (unsigned)'!');
}

/* produce pileup data (calls, quals, and read depths) from current
   position for sample s. manage buffer reallocation of pd fields. */
void
pileup_current_data(unsigned s, struct pileup_data *pd)
{
    unsigned n_base_ct;
    struct bqs_count *base_ct = NULL;
    pileup_current_bqs(s, &base_ct, &n_base_ct);
    unsigned b, n_calls = 0;
    pd->n_match_lo_q = 0;
    pd->n_match_hi_q = 0;

    for (b = 0; b != n_base_ct; ++b)
        if (base_ct[b].qual < min_quality_score)
            pd->n_match_lo_q += base_ct[b].ct;
        else
            pd->n_match_hi_q += base_ct[b].ct;

    n_calls = pd->n_match_lo_q + pd->n_match_hi_q;
    pd->quals.size = n_calls;
    ALLOC_GROW(pd->quals.buf, pd->quals.size, pd->quals.alloc);

    struct indel_count *indel_ct = NULL;
    unsigned n_indel_ct, indel_len_total = 0;

    pileup_current_indels(s, &indel_ct, &n_indel_ct);
    struct indel_seq *isq;
    pd->n_indel = 0;
    for (b = 0; b != n_indel_ct; ++b) {
        isq = kh_key(tls.indel_hash, indel_ct[b].indel_itr);
        indel_len_total += indel_ct[b].ct * (strlen(isq->seq) + 10); /* reasonable upper bound */
        pd->n_indel += indel_ct[b].ct;
    }
    pd->calls.size = n_calls + indel_len_total;
    ALLOC_GROW(pd->calls.buf, pd->calls.size, pd->calls.alloc);
    
    /* print out the calls */
    struct pileup_locus_info pli;
    pileup_current_info(&pli);
    unsigned refbase_i = get_cur_refbase_code16();
    unsigned p = 0, p_cur, p_end;
    for (b = 0; b != n_base_ct; ++b) {
        p_end = p + base_ct[b].ct;
        char bc = bqs_count_to_call(base_ct[b], refbase_i);
        char qs = qual_to_char(base_ct[b].qual);
        p_cur = p;
        for (; p != p_end; ++p) pd->calls.buf[p] = bc;
        for (p = p_cur; p != p_end; ++p) pd->quals.buf[p] = qs;
    }
    pd->quals.size = p;

    /* print out all indel representations */
    for (b = 0; b != n_indel_ct; ++b) {
        isq = kh_key(tls.indel_hash, indel_ct[b].indel_itr);
        pd->calls.buf[p++] = isq->is_ins ? '+' : '-';
        p += sprintf(pd->calls.buf + p, "%Zu", strlen(isq->seq));
        strcpy(pd->calls.buf + p, isq->seq);
        p += strlen(isq->seq);
    }
    pd->calls.size = p;
}


void
free_pileup_data(struct pileup_data *pd)
{
    if (pd->calls.buf != NULL) {
        free(pd->calls.buf);
        pd->calls.alloc = 0;
    }
    if (pd->quals.buf != NULL) {
        free(pd->quals.buf);
        pd->quals.alloc = 0;
    }
}


/* return the complement of the 4-bit integer representation of the
   base. this happens to be the bits in reverse order. */
static inline int
cmpl(int base)
{
    static int rev[] = {
        0, 15, 4, 12, 2, 10, 6, 14, 1, 9, 5, 13, 3, 11, 7, 15
    };
    return rev[base];
}


/* traverse b, tallying the match blocks into ts->pbqt_hash, and the
   indels into ts->indel_ct_hash and tls.indel_hash as necessary */
static void
process_bam_stats(bam1_t *b, struct tally_stats *ts)
{
    int32_t tid = b->core.tid, qpos = 0, rpos = b->core.pos, q, r;
    unsigned c;
    uint32_t *cigar = bam_get_cigar(b), op, ln;
    unsigned strand = bam_is_rev(b) ? 0 : 1;
    int ret;
    khash_t(pbqt_h) *ph = ts->pbqt_hash;
    uint8_t *bam_seq = bam_get_seq(b);
    uint8_t *bam_qual = bam_get_qual(b);
    
    for (c = 0; c != b->core.n_cigar; ++c) {
        op = bam_cigar_op(cigar[c]);
        ln = bam_cigar_oplen(cigar[c]);

        khiter_t it;
        int32_t qe;
        if (op == BAM_CMATCH || op == BAM_CEQUAL || op == BAM_CDIFF) {
            if (0) {
                int32_t qb = b->core.l_qseq - qpos - 1;
                qe = qb - ln;
                for (q = qb, r = rpos; q != qe; --q, ++r) {
                    union pbqt_key stat = {
                        .v = { r, tid, cmpl(bam_seqi(bam_seq, q)), bam_qual[q], strand }
                    };
                    if ((it = kh_get(pbqt_h, ph, stat.k)) == kh_end(ph)) {
                        it = kh_put(pbqt_h, ph, stat.k, &ret);
                        kh_val(ph, it) = 0;
                    }
                    kh_val(ph, it)++;
                }
            } else {
                qe = qpos + ln;
                for (q = qpos, r = rpos; q != qe; ++q, ++r) {
                    union pbqt_key stat = {
                        .v = { r, tid, bam_seqi(bam_seq, q), bam_qual[q], strand }
                    };
                    if ((it = kh_get(pbqt_h, ph, stat.k)) == kh_end(ph)) {
                        it = kh_put(pbqt_h, ph, stat.k, &ret);
                        kh_val(ph, it) = 0;
                    }
                    kh_val(ph, it)++;
                }
            }
        } else if (op == BAM_CINS || op == BAM_CDEL) {
            /* tally insertion of query. */
            struct indel_seq *isq = malloc(sizeof(struct indel_seq) + ln + 1);
            isq->is_ins = (op == BAM_CINS);
            unsigned i;

            if (isq->is_ins) {
                if (bam_is_rev(b)) {
                    int32_t qb = b->core.l_qseq - qpos - 1;
                    for (i = 0, q = qb; i != ln; ++i, --q)
                        isq->seq[i] = seq_nt16_str[cmpl(bam_seqi(bam_seq, q))];
                } else {
                    for (i = 0, q = qpos; i != ln; ++i, ++q)
                        isq->seq[i] = seq_nt16_str[bam_seqi(bam_seq, q)];
                }
            }
            else
                memcpy(isq->seq, reference_seq.contig[tid].seq + rpos, ln);

            isq->seq[ln] = '\0';

            khiter_t it;
            if ((it = kh_get(indel_h, tls.indel_hash, isq)) == kh_end(tls.indel_hash))
                it = kh_put(indel_h, tls.indel_hash, isq, &ret);
            else free(isq);
            
            struct contig_pos cpos = { tid, rpos };
            incr_indel_count(ts, cpos, it);
        } else
            ; /* All other operations (N, S, H, P) do not result in
                 tallying anything */
        
        /* but, we do advance the query */
        if (bam_cigar_type(op) & 1) qpos += ln;
        if (bam_cigar_type(op) & 2) rpos += ln;
    }
}


/* process an entire set of raw (uncompressed) BAM records in [rec,
   end).  maintains a map of possibly overlapping reads and resolves
   their quality scores. returns the last position processed, or { 0,
   0 } if given empty input. */
struct contig_pos
process_bam_block(char *rec, char *end, struct tally_stats *ts)
{
    bam1_t b, b_mate;
    int ret;
    char *rec_next;
    struct contig_pos last_pos = { 0, 0 };
    while (rec != end) {
        rec_next = bam_parse(rec, &b);
        if (! (b.core.flag & BAM_FPAIRED) || ! (b.core.flag & BAM_FPROPER_PAIR)) {
            /* either this read is not paired or it is paired but not
               mapped in a proper pair.  tally and don't store. */
            process_bam_stats(&b, ts);
        } else { /* b is paired and mapped in proper pair */
            if (b.core.pos < b.core.mpos) { /* b is upstream mate */
                if (bam_endpos(&b) < b.core.mpos) { /* does not overlap with mate */
                    process_bam_stats(&b, ts);
                } else {
                    /* b overlaps mate.  since b is upstream, it won't
                       have been stored in overlaps hash yet. store;
                       do not tally yet. */
                    khiter_t it =
                        kh_put(olap_h, ts->overlap_hash, bam_get_qname(&b), &ret);
                    kh_val(ts->overlap_hash, it) = rec;
                }
            } else {
                /* b is downstream mate. must search overlaps hash to
                   see if it overlaps with its upstream mate. */
                khiter_t it =
                    kh_get(olap_h, ts->overlap_hash, bam_get_qname(&b));
                if (it != kh_end(ts->overlap_hash)) {
                    (void)bam_parse(kh_val(ts->overlap_hash, it), &b_mate);
                    // !!! NEED TO IMPLEMENT: tweak_overlap_quality(&b, &b_mate);
                    process_bam_stats(&b, ts);
                    process_bam_stats(&b_mate, ts);
                    kh_del(olap_h, ts->overlap_hash, it);
                } else {
                    if (b.core.pos != b.core.mpos) {
                        /* downstream mate does not overlap with its
                           upstream partner */
                        process_bam_stats(&b, ts);
                    } else {
                        /* b is first in the pair to be encountered.  store it */
                        khiter_t it = 
                            kh_put(olap_h, ts->overlap_hash, bam_get_qname(&b), &ret);
                        kh_val(ts->overlap_hash, it) = rec;
                    }
                }
            }
        }
        if (rec_next == end) last_pos = (struct contig_pos){ b.core.tid, b.core.pos };
        rec = rec_next;
    }
    return last_pos;
}


/* marginalize out q and t from pbqt_hash, storing results in pb_hash.
   counts are only tallied if q >= min_quality_score (global var).   */
static khash_t(pb_h) *
summarize_base_counts(unsigned s, struct contig_pos tally_end)
{
    union pbqt_key src_k;
    khash_t(pbqt_h) *src_h = tls.ts[s].pbqt_hash;

    union pos_key trg_k;
    khash_t(pb_h) *trg_h = kh_init(pb_h);

    unsigned ct;
    kh_clear(pb_h, trg_h);
    khiter_t src_itr, trg_itr;
    int ret;
    for (src_itr = kh_begin(src_h); src_itr != kh_end(src_h); ++src_itr) {
        if (! kh_exist(src_h, src_itr)) continue;
        src_k.k = kh_key(src_h, src_itr);
        trg_k.v = (struct contig_pos){ src_k.v.tid, src_k.v.pos };
        if (! less_contig_pos(trg_k.v, tally_end)) continue;

        ct = kh_val(src_h, src_itr);
        trg_itr = kh_get(pb_h, trg_h, trg_k.k);
        if (trg_itr == kh_end(trg_h)) {
            trg_itr = kh_put(pb_h, trg_h, trg_k.k, &ret);
            assert(ret != 0);
            kh_val(trg_h, trg_itr) = 
                (struct base_count){ .ct_filt = { 0, 0, 0, 0 }, 
                                     .n_match_lo_q = 0,
                                     .n_match_hi_q = 0,
                                     .n_match_fuzzy = 0 };
        }
        int pure_base = seq_nt16_int[src_k.v.base];
        if (pure_base == 4) 
            kh_val(trg_h, trg_itr).n_match_fuzzy += ct;
        else {
            if (src_k.v.qual >= min_quality_score) {
                kh_val(trg_h, trg_itr).ct_filt[pure_base] += ct;
                kh_val(trg_h, trg_itr).n_match_hi_q += ct;
            }
            else
                kh_val(trg_h, trg_itr).n_match_lo_q += ct;
        }
    }
    return trg_h;
}


static int
pos_base_count_less(struct pos_base_count a, struct pos_base_count b)
{
    return less_contig_pos(a.cpos, b.cpos);
}

KSORT_INIT(pbc_sort, struct pos_base_count, pos_base_count_less);

/* condense base counts information in pbqt hash into a sorted
   array */
void
pileup_prepare_basecalls(unsigned s)
{
    khash_t(pb_h) *ph = summarize_base_counts(s, tls.tally_end);
    unsigned n = kh_size(ph);
    struct tally_stats *ts = &tls.ts[s];
    struct pos_base_count *ary = 
        realloc(ts->base_ct, n * sizeof(struct pos_base_count));

    khiter_t it;
    unsigned i;
    union pos_key k;
    for (it = kh_begin(ph), i = 0; it != kh_end(ph); ++it)
        if (kh_exist(ph, it)) {
            k.k = kh_key(ph, it);
            ary[i++] = (struct pos_base_count){ k.v, kh_val(ph, it) };
        }
    
    ks_introsort(pbc_sort, n, ary);
    ts->base_ct = ary;
    ts->base_cur = ary;
    ts->base_end = ary + i;
}

static int
pos_iter_indel_count_less(struct pos_indel_count a, struct pos_indel_count b)
{
    return less_contig_pos(a.cpos, b.cpos)
        || (equal_contig_pos(a.cpos, b.cpos)
            && a.ict.indel_itr < b.ict.indel_itr);
}

KSORT_INIT(pi_sort, struct pos_indel_count, pos_iter_indel_count_less);


/* extract indel data in indel_ct_hash into an array sorted by
   position, for fast retrieval of the current position's indel
   data. */
void
pileup_prepare_indels(unsigned s)
{
    struct tally_stats *ts = &tls.ts[s];
    khash_t(indel_ct_h) *ih = ts->indel_ct_hash;
    unsigned n = ts->n_indel_nodes;
    struct pos_indel_count *ary = 
        realloc(ts->indel_ct, n * sizeof(struct pos_indel_count));
    
    khiter_t it;
    unsigned i;
    struct indel_count_node *node;
    union pos_key k;
    for (it = kh_begin(ih), i = 0; it != kh_end(ih); ++it)
        if (kh_exist(ih, it)) {
            k.k = kh_key(ih, it);
            if (less_contig_pos(k.v, tls.tally_end)) {
                node = kh_val(ih, it);
                while (node) {
                    ary[i++] = (struct pos_indel_count){ k.v, node->ict };
                    node = node->next;
                }
            }
        }

    ks_introsort(pi_sort, n, ary);
    ts->indel_ct = ary;
    ts->indel_cur = ary;
    ts->indel_end = ary + i;
}

/* */
static void
free_indel_counts(struct indel_count_node *nd)
{
    struct indel_count_node *t;
    while (nd) {
        t = nd;
        nd = nd->next;
        free(t);
    }
}


/* increment a count of a particular indel, represented by indel_key,
   in a sample s at pos. */
static void
incr_indel_count(struct tally_stats *ts,
                 struct contig_pos pos,
                 khiter_t indel_itr)
{
    union pos_key k = { .v = pos };
    khash_t(indel_ct_h) *ih = ts->indel_ct_hash;

    khiter_t i;
    int ret;
    if ((i = kh_get(indel_ct_h, ih, k.k)) == kh_end(ih)) {
        i = kh_put(indel_ct_h, ih, k.k, &ret);
        kh_val(ih, i) = NULL;
    }

    struct indel_count_node *head, *node;
    head = node = kh_val(ih, i);
    
    /* find node with matching key, or NULL */
    while (node && node->ict.indel_itr != indel_itr)
        node = node->next;

    if (node) ++node->ict.ct;
    else {
        ++ts->n_indel_nodes;
        node = malloc(sizeof(struct indel_count_node));
        *node = (struct indel_count_node){ head, { indel_itr, 1 } };
        kh_val(ih, i) = node;
    }
}




/* clear statistics that are no longer needed. call this function
   after pileup_next_pos() returns 1. */
void
pileup_clear_finished_stats()
{
    khiter_t it;
    unsigned s;
    struct contig_pos kpos;

    for (s = 0; s != tls.n_samples; ++s) {
        /* clear finished entries in pbqt hash */
        struct tally_stats *ts = &tls.ts[s];
        khash_t(pbqt_h) *ph = ts->pbqt_hash;
        union pbqt_key pk;
        for (it = kh_begin(ph); it != kh_end(ph); ++it) {
            if (! kh_exist(ph, it)) continue;
            pk.k = kh_key(ph, it);
            kpos = (struct contig_pos){ pk.v.tid, pk.v.pos };
            if (less_contig_pos(kpos, tls.tally_end))
                kh_del(pbqt_h, ph, it);
        }

        /* clear finished entries in indel counts hash */
        khash_t(indel_ct_h) *ih = ts->indel_ct_hash;
        union pos_key pk2;
        for (it = kh_begin(ih); it != kh_end(ih); ++it) {
            if (! kh_exist(ih, it)) continue;
            pk2.k = kh_key(ih, it);
            if (less_contig_pos(pk2.v, tls.tally_end)) {
                free_indel_counts(kh_val(ih, it));
                kh_del(indel_ct_h, ih, it);
            }
        }

        /* completely clear temporary data */
        free(ts->base_ct);
        ts->base_ct = NULL;
        free(ts->bqs_ct);
        ts->bqs_ct = NULL;
        free(ts->indel_ct);
        ts->indel_ct = NULL;
    }
}


/* advance tls.cur_pos to next position for which at least one sample
   has base or indel entries. update base_cur and indel_cur pointers
   for each sample so that calls to pileup_{basecall,bqs,indel}_stats
   return the right values. return 1 if there is a next position, 0 if
   reached the end. if this is the first call of all time, or the
   first call after a previous end, return the first available
   position. */
int
pileup_next_pos()
{
    unsigned s;
    struct tally_stats *ts;
    
    /* update all cur pointers.  */
    struct contig_pos suc;
    if (tls.cur_pos.tid == TID_UNSET)
        suc = (struct contig_pos){ 0, 0 };
    else
        suc = (struct contig_pos){ tls.cur_pos.tid, tls.cur_pos.pos + 1 };

    for (s = 0; s != tls.n_samples; ++s) {
        ts = &tls.ts[s];
        /* increment bqs iterator ('while' is used here, because bqs
           have multiple entries per position) */
        while (ts->bqs_cur != ts->bqs_end
            && less_contig_pos(ts->bqs_cur->cpos, suc))
            ++ts->bqs_cur;

        /* increment base iterator */
        if (ts->base_cur != ts->base_end
            && less_contig_pos(ts->base_cur->cpos, suc))
            ++ts->base_cur;

        /* increment indel iterator ('while' used here, because may be
           multiple entries per position) */
        while (ts->indel_cur != ts->indel_end
            && less_contig_pos(ts->indel_cur->cpos, suc))
            ++ts->indel_cur;
    }
    /* recalculate tls.cur_pos to be the minimum position among all
       valid base_cur and indel_cur pointers in all samples. (we don't
       need to test ts->bqs_cur since it is */
    struct contig_pos unset_pos = { TID_UNSET, 0 };
    struct contig_pos tmp_b, tmp_i, tmp_q, tmp = unset_pos;
    for (s = 0; s != tls.n_samples; ++s) {
        ts = &tls.ts[s];
        tmp_q = (ts->bqs_cur == ts->bqs_end) ? tmp : ts->bqs_cur->cpos;
        tmp_b = (ts->base_cur == ts->base_end) ? tmp : ts->base_cur->cpos;
        tmp_i = (ts->indel_cur == ts->indel_end) ? tmp : ts->indel_cur->cpos;
        tmp = MIN_CONTIG_POS(tmp, tmp_q);
        tmp = MIN_CONTIG_POS(tmp, tmp_b);
        tmp = MIN_CONTIG_POS(tmp, tmp_i);
    }
    tls.cur_pos = tmp;
    return less_contig_pos(tls.cur_pos, unset_pos);
}


void
pileup_current_info(struct pileup_locus_info *pli)
{
    unsigned tid = tls.cur_pos.tid;
    strcpy(pli->refname, reference_seq.contig[tid].name);
    pli->refbase = reference_seq.contig[tid].seq[tls.cur_pos.pos];
    pli->pos = tls.cur_pos.pos;
}


void
pileup_get_indel(unsigned indel_id, struct indel_seq **indel)
{
    khiter_t k = indel_id;
    struct indel_seq *is = kh_key(tls.indel_hash, k);
    size_t sz = sizeof(struct indel_seq) + strlen(is->seq);
    *indel = realloc(*indel, sz);
    memcpy(*indel, is, sz);
}


void
pileup_final_input()
{
    tls.tally_end = (struct contig_pos){ UINT_MAX, UINT_MAX };
}
