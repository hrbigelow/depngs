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
    uint32_t tid: 22;
    uint32_t base: 2;
    uint32_t qual: 7;
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
    unsigned no_more_input; /* set to 1 to signal */
    khash_t(olap_h) *overlap_hash;

    /* summary statistics are compiled for positions < tally_end */
    khash_t(pb_h) *pb_hash;
    struct pos_base_count *base_ct, *base_cur, *base_end;
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



/* */
void
batch_pileup_thread_init(unsigned n_samples)
{
    tls.indel_hash = kh_init(indel_h);
    tls.n_samples = n_samples;
    tls.tally_end = (struct contig_pos){ 0, 0 };
    tls.ts = malloc(n_samples * sizeof(struct tally_stats));
    unsigned s;
    for (s = 0; s != n_samples; ++s)
        tls.ts[s] = (struct tally_stats){
            .pbqt_hash = kh_init(pbqt_h),
            .indel_ct_hash = kh_init(indel_ct_h),
            .n_indel_nodes = 0,
            .no_more_input = 0,
            .overlap_hash = kh_init(olap_h),
            .pb_hash = kh_init(pb_h),
            .base_ct = NULL, 
            .base_cur = NULL, 
            .base_end = NULL,
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
            free(kh_val(tls.indel_hash, it));

    kh_destroy(indel_h, tls.indel_hash);

    unsigned s;
    for (s = 0; s != tls.n_samples; ++s) {
        struct tally_stats *ts = &tls.ts[s];
        kh_destroy(pbqt_h, ts->pbqt_hash);
        kh_destroy(indel_ct_h, ts->indel_ct_hash);
        kh_destroy(olap_h, ts->overlap_hash);
        kh_destroy(pb_h, ts->pb_hash);
        free(ts->base_ct);
        free(ts->indel_ct);
    }
}

static void
summarize_base_counts(unsigned s, struct contig_pos tally_end);

static void
make_p_array(unsigned s);

static void
make_indel_array(unsigned s, struct contig_pos tally_end);

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
   sample. update tls.tally_end to reflect furthest position seen. */
void
tally_pileup_stats(const struct managed_buf bam, unsigned s)
{
    char *b = bam.buf, *be = b + bam.size;
    struct contig_pos bpos = 
        process_bam_block(b, be, &tls.ts[s]);
    tls.tally_end = MAX_CONTIG_POS(tls.tally_end, bpos);
}


/* summarize statistics for sample s, up to tls.tally_end, or to
   completion if there is no more input available */
void
summarize_pileup_stats(unsigned s)
{
    static struct contig_pos last = { UINT_MAX, UINT_MAX };
    struct contig_pos end = 
        tls.ts[s].no_more_input 
        ? last
        : tls.tally_end;
    
    summarize_base_counts(s, end);
    make_p_array(s);
    make_indel_array(s, end);
}


/* return 0,1,2,3 if the current position reference base is A,C,G,T
   (or their lowercase equivalent), or 4 for anything else. */
static unsigned
get_cur_refbase_index()
{
    char refbase = reference_seq.contig[tls.cur_pos.tid].seq[tls.cur_pos.pos];
    static char nucs[] = "ACGT";
    char *p;
    return (p = index(nucs, toupper(refbase))) == NULL ? 4 : p - nucs;
}


/* provide basecall stats for a sample at current position, or the
   null statistic. */
struct base_count
pileup_basecall_stats(unsigned s)
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
        return refsam_ct[get_cur_refbase_index()];

    struct tally_stats *ts = &tls.ts[s];
    if (ts->base_cur == ts->base_end
        || less_contig_pos(tls.cur_pos, ts->base_cur->cpos))
        return null_ct;

    else
        return ts->base_cur->bct;
}


#define MIN_QUAL 0
#define MAX_QUAL 127

/* provide (b,q,s) stats for a sample at current position.  *cts is
   reallocated as necessary. *n_cts set to number of distinct stats
   that are populated. complete statistics are provided; there is no
   filtering by quality score. */
void
pileup_bqs_stats(unsigned s, struct bqs_count **cts, unsigned *n_cts)
{
    unsigned b, q, st;
    union pbqt_key pk;
    khint_t j;
    pk.v.tid = tls.cur_pos.tid;
    pk.v.pos = tls.cur_pos.pos;
    unsigned n = 0;

    static struct bqs_count ref_bqs[] = {
        { 0, 50, 0, 1e6 },
        { 1, 50, 0, 1e6 },
        { 2, 50, 0, 1e6 },
        { 3, 50, 0, 1e6 }
    };

    if (s == REFERENCE_SAMPLE) {
        *n_cts = 1;
        *cts = realloc(*cts, 1 * sizeof(struct bqs_count));
        (*cts)[0] = ref_bqs[get_cur_refbase_index()];
        return;
    }
        
    /* should be about 4k.  is this too large? */
    struct bqs_count buf[(MAX_QUAL - MIN_QUAL) * 4 * 2];
    
    struct tally_stats *ts = &tls.ts[s];
    for (st = 0; st != 2; ++st) {
        pk.v.strand = st;
        for (b = 0; b != 4; ++b) {
            pk.v.base = b;
            for (q = MIN_QUAL; q != MAX_QUAL; ++q) {
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
void
pileup_indel_stats(unsigned s, struct indel_count **cts, unsigned *n_cts)
{
    if (s == REFERENCE_SAMPLE) {
        *n_cts = 0;
        return;
    }
    
    /* pre-scan to find total number of distinct indels (use linear
       scan since this will be very small) */
    struct pos_indel_count *pc = tls.ts[s].indel_cur;
    unsigned n = 0;
    while (pc != tls.ts[s].indel_end &&
           equal_contig_pos(pc->cpos, tls.cur_pos))
        ++n;

    *cts = realloc(*cts, n * sizeof(struct indel_count));
    *n_cts = n;
    unsigned i = 0;
    pc = tls.ts[s].indel_cur;
    while (i != n)
        *cts[i++] = pc++->ict;
}


/* produce pileup strings from current position for sample s,
   storing in call and qual, respectively */
void
pileup_pair_indel_stats(struct indel_count *cts1, unsigned n_cts1,
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


static inline char
bqs_count_to_call(struct bqs_count bc, unsigned refbase_i)
{
    static char *calls[] = { "ACGT.", "acgt," };
    return calls[bc.strand][bc.base == refbase_i ? 4 : bc.base];
}

/* convert a BAM quality to a character */
static inline char
qual_to_char(unsigned q)
{
    return (char)(q - 33);
}

/* produce pileup data (calls, quals, and read depths) from current
   position for sample s */
void
pileup_current_data(unsigned s, struct pileup_data *pd)
{
    unsigned n_base_ct;
    struct bqs_count *base_ct = NULL, *bcp;
    pileup_bqs_stats(s, &base_ct, &n_base_ct);
    unsigned b, n_calls = 0;
    pd->n_match_lo_q = 0;
    pd->n_match_hi_q = 0;

    for (b = 0; b != n_base_ct; ++b) {
        bcp = base_ct + b;
        if (bcp->qual < min_quality_score)
            pd->n_match_lo_q += bcp->ct;
        else
            pd->n_match_hi_q += bcp->ct;
    }

    pd->quals.size = pd->n_match_lo_q + pd->n_match_hi_q;
    ALLOC_GROW(pd->quals.buf, pd->quals.size, pd->quals.alloc);

    struct indel_count *indel_ct = NULL;
    unsigned n_indel_ct, indel_len_total = 0;

    pileup_indel_stats(s, &indel_ct, &n_indel_ct);
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
    static char nucs[] = "ACGT";
    unsigned refbase_i = index(nucs, pli.refbase) - nucs;
    unsigned p = 0, p_end;
    for (b = 0; b != n_base_ct; ++b) {
        p_end = p + base_ct[b].ct;
        char bc = bqs_count_to_call(base_ct[b], refbase_i);
        char qs = qual_to_char(base_ct[b].qual);
        for (; p != p_end; ++p) pd->calls.buf[p] = bc;
        for (; p != p_end; ++p) pd->quals.buf[p] = qs;
    }
    pd->quals.buf[p] = '\0';
    pd->quals.size = p;

    /* print out all indel representations */
    for (b = 0; b != n_indel_ct; ++b) {
        isq = kh_key(tls.indel_hash, indel_ct[b].indel_itr);
        pd->calls.buf[p++] = isq->is_ins ? '+' : '-';
        p += sprintf(pd->calls.buf + p, "%Zu", strlen(isq->seq));
        strcpy(pd->calls.buf + p, isq->seq);
        p += strlen(isq->seq);
    }
    pd->calls.buf[p] = '\0';
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
    
    for (c = 0; c != b->core.n_cigar; ++c) {
        op = bam_cigar_op(cigar[c]);
        ln = bam_cigar_oplen(cigar[c]);

        if (op == BAM_CMATCH || op == BAM_CEQUAL || op == BAM_CDIFF) {
            int32_t qe = qpos + ln;
            khiter_t it;
            for (q = qpos, r = rpos; q != qe; ++q, ++r) {
                union pbqt_key stat = {
                    .v = {
                        r, tid, 
                        bam_seqi(bam_get_seq(b), q), 
                        bam_get_qual(b)[q],
                        strand 
                    }
                };
                if ((it = kh_get(pbqt_h, ph, stat.k)) == kh_end(ph)) {
                    it = kh_put(pbqt_h, ph, stat.k, &ret);
                    kh_val(ph, it) = 0;
                }
                kh_val(ph, it)++;
            }
        } else if (op == BAM_CINS || op == BAM_CDEL) {
            /* tally insertion of query. */
            struct indel_seq *isq = malloc(sizeof(struct indel_seq) + ln + 1);
            isq->is_ins = (op == BAM_CINS);
            unsigned i;

            if (isq->is_ins)
                for (i = 0, q = qpos; i != ln; ++i, ++q)
                    isq->seq[i] = seq_nt16_str[bam_seqi(bam_get_seq(b), q)];

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
static void
summarize_base_counts(unsigned s, struct contig_pos tally_end)
{
    union pbqt_key src_k;
    khash_t(pbqt_h) *src_h = tls.ts[s].pbqt_hash;

    union pos_key trg_k;
    khash_t(pb_h) *trg_h = tls.ts[s].pb_hash;

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
            kh_val(trg_h, trg_itr) = (struct base_count){ { 0, 0, 0, 0 }, 0, 0 };
        }
        if (src_k.v.qual >= min_quality_score) {
            kh_val(trg_h, trg_itr).ct_filt[src_k.v.base] += ct;
            kh_val(trg_h, trg_itr).n_match_hi_q += ct;
        } else {
            kh_val(trg_h, trg_itr).n_match_lo_q += ct;
        }
    }
}


static int
pos_base_count_less(struct pos_base_count a, struct pos_base_count b)
{
    return less_contig_pos(a.cpos, b.cpos);
}

KSORT_INIT(pbc_sort, struct pos_base_count, pos_base_count_less);

/* create a sorted array from sample s's pb_hash.  Note: pb_hash is only
   summarized for all positions < tally_end). also initialize base_ct,
   base_cur, and base_end for this sample. */
static void
make_p_array(unsigned s)
{
    khash_t(pb_h) *ph = tls.ts[s].pb_hash;
    unsigned n = kh_size(ph);
    struct pos_base_count *ary = 
        realloc(tls.ts[s].base_ct, n * sizeof(struct pos_base_count));

    khiter_t it;
    unsigned i;
    union pos_key k;
    for (it = kh_begin(ph), i = 0; it != kh_end(ph); ++it)
        if (kh_exist(ph, it)) {
            k.k = kh_key(ph, it);
            ary[i++] = (struct pos_base_count){ k.v, kh_val(ph, it) };
        }
    
    ks_introsort(pbc_sort, n, ary);
    tls.ts[s].base_ct = ary;
    tls.ts[s].base_cur = ary;
    tls.ts[s].base_end = ary + i;
}

static int
pos_iter_indel_count_less(struct pos_indel_count a, struct pos_indel_count b)
{
    return less_contig_pos(a.cpos, b.cpos)
        || (equal_contig_pos(a.cpos, b.cpos)
            && a.ict.indel_itr < b.ict.indel_itr);
}

KSORT_INIT(pi_sort, struct pos_indel_count, pos_iter_indel_count_less);


/* create a sorted array from sample s's indel_hash */
static void
make_indel_array(unsigned s, struct contig_pos tally_end)
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
        if (kh_exist(ih, it)) {
            k.k = kh_key(ih, it);
            if (less_contig_pos(k.v, tally_end)) {
                node = kh_val(ih, it);
                while (node) {
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
static void
free_indel_counts(struct indel_count_node *nd)
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
   after next_base() and next_indel() return 1. */
void
pileup_clear_finished_stats()
{
    khiter_t it;
    unsigned s;
    struct contig_pos kpos;

    for (s = 0; s != tls.n_samples; ++s) {
        /* clear finished entries in pbqt hash */
        khash_t(pbqt_h) *ph = tls.ts[s].pbqt_hash;
        union pbqt_key pk;
        for (it = kh_begin(ph); it != kh_end(ph); ++it) {
            if (! kh_exist(ph, it)) continue;
            pk.k = kh_key(ph, it);
            kpos = (struct contig_pos){ pk.v.tid, pk.v.pos };
            if (less_contig_pos(kpos, tls.tally_end))
                kh_del(pbqt_h, ph, it);
        }

        /* clear finished entries in indel counts hash */
        khash_t(indel_ct_h) *ih = tls.ts[s].indel_ct_hash;
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
        kh_clear(pb_h, tls.ts[s].pb_hash);
        free(tls.ts[s].base_ct);
        free(tls.ts[s].indel_ct);
    }
}


/* advance tls.cur_pos to next position for which at least one sample
   has base or indel entries. update base_cur and indel_cur pointers
   for each sample so that calls to pileup_{basecall,bqs,indel}_stats
   return the right values. return 1 if there is a next position, 0 if
   reached the end. */
int
pileup_next_pos()
{
    unsigned s;
    struct tally_stats *ts;
    
    /* update all cur pointers */
    struct contig_pos suc =
        { tls.cur_pos.tid, tls.cur_pos.pos + 1 };
    for (s = 0; s != tls.n_samples; ++s) {
        ts = &tls.ts[s];
        /* increment base iterator */
        if (ts->base_cur != ts->base_end
            && less_contig_pos(ts->base_cur->cpos, suc))
            ++ts->base_cur;
        
        /* increment indel iterator */
        if (ts->indel_cur != ts->indel_end
            && less_contig_pos(ts->indel_cur->cpos, suc))
            ++ts->indel_cur;
        
    }
    /* recalculate tls.cur_pos to be the minimum position of either an  */
    struct contig_pos end_pos = { UINT_MAX, UINT_MAX }, 
        tmp = end_pos, tmp_b, tmp_i;
        for (s = 0; s != tls.n_samples; ++s) {
            ts = &tls.ts[s];
            tmp_b = (ts->base_cur == ts->base_end) ? tmp : ts->base_cur->cpos;
            tmp_i = (ts->indel_cur == ts->indel_end) ? tmp : ts->indel_cur->cpos;
            tmp = MIN_CONTIG_POS(tmp, tmp_b);
            tmp = MIN_CONTIG_POS(tmp, tmp_i);
        }
        tls.cur_pos = tmp;
        return less_contig_pos(tls.cur_pos, end_pos);
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
