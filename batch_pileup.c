/* 
   Populates complete pileup statistics for multiple samples in
   tandem.  Provides access to these statistics for a given sample at
   the current position in three views: basecall counts, (b,q,s)
   triplet counts, or indel counts.  Maintains an internal 'iterator'
   position and provides the user with a 'next' functionality for the
   iterator.
*/

#include "batch_pileup.h"
#include "ksort.h"
#include "cache.h"
#include "bam_reader.h"
#include "bam_sample_info.h"
#include "fasta.h"
#include "locus_range.h"

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

union pos_key {
    khint64_t k;
    struct contig_pos v;
};


struct pos_base_count {
    struct contig_pos cpos;
    struct base_count bct;
};

struct pos_bqs_count {
    struct contig_pos cpos;
    struct bqs_count bqs_ct;
};

struct pos_indel_count {
    struct contig_pos cpos;
    struct indel_count ict;
};


static int
less_indel(struct indel a, struct indel b)
{
    return a.is_ins < b.is_ins
        || (a.is_ins == b.is_ins
            && a.length < b.length);
}



/* hash type for storing bqt (basecall, quality, strand) counts at
   position.  uses 'union pbqt_key' as the key  */
KHASH_MAP_INIT_INT64(pbqt_h, unsigned);

/* hash type for storing basecall counts at position. uses 'union
   pos_key' as the key. */
KHASH_MAP_INIT_INT64(pb_h, struct base_count);

/* hash type for storing plain old strings (for the insertion
   sequences) */
KHASH_SET_INIT_STR(str_h);

/* hash type for tallying indel events. use   */
KHASH_MAP_INIT_INT64(indel_ct_h, struct indel_count_ary);

/* overlapping mates hash, key is q_name.  values are strdup-allocated
   raw BAM records. */
KHASH_MAP_INIT_STR(olap_h, char *);


/* complete tally statistics for the set of BAM records overlapping a
   given set of locus ranges. */
struct tally_stats {
    /* Once tallying is complete, pbqt_hash and indel_ct_hash have
       complete statistics for loci < tally_end, and
       partial statistics for loci in [tally_end, MAX). */
    khash_t(pbqt_h) *pbqt_hash;
    khash_t(indel_ct_h) *indel_ct_hash;
    khash_t(olap_h) *overlap_hash;

    /* summary statistics are compiled for positions < tally_end */
    struct pos_base_count *base_ct, *base_cur, *base_end;
    struct pos_bqs_count *bqs_ct, *bqs_cur, *bqs_end;
    struct pos_indel_count *indel_ct, *indel_cur, *indel_end;
};



static __thread struct tls {
    khash_t(str_h) *seq_hash;
    unsigned n_samples;
    struct contig_pos tally_end;
    struct tally_stats *ts; /* tally stats for each sample */
    struct contig_pos cur_pos; /* current position being queried */

    /* these fields refreshed at pileup_tally_stats, freed at
       pileup_clear_finished_stats().  they define the region of
       interest. */
    struct contig_fragment *refseqs, *cur_refseq;
    unsigned n_refseqs;
    struct base_count null_ct;
    struct base_count refsam_ct[5];
} tls;


#define FRAGMENT_MAX_INLINE_SEQLEN 10
struct contig_fragment {
    struct contig_region reg;
    char buf[FRAGMENT_MAX_INLINE_SEQLEN + 1];
    char *seq;
};


static unsigned pseudo_depth;
static unsigned min_quality_score;
static unsigned skip_empty_loci; /* if 1, pileup_next_pos() advances
                                      past loci that have no data. */


void
batch_pileup_init(unsigned min_qual, unsigned skip_empty,
                  unsigned _pseudo_depth)
{
    min_quality_score = min_qual;
    skip_empty_loci = skip_empty;
    pseudo_depth = _pseudo_depth;
}


void
batch_pileup_free()
{
}

/* used to represent the position after we run out of loci in a
   chunk. also defined as the maximum possible position. tls.cur_pos
   is set to this value when pileup_next_pos() is called and
   tls.cur_pos is incremented out of the region of interest */
static const struct contig_pos g_end_pos = { UINT_MAX, UINT_MAX };

/* use to represent the current position after a new batch was loaded
   but before pileup_next_pos() is called. tls.cur_pos is set to this
   value after calling pileup_clear_stats() and also after
   batch_pileup_thread_init() */
static const struct contig_pos g_unset_pos = { UINT_MAX - 1, 0 };


/* */
void
batch_pileup_thread_init(unsigned n_samples, const char *fasta_file)
{
    fasta_thread_init(fasta_file);

    tls.seq_hash = kh_init(str_h);
    tls.n_samples = n_samples;
    tls.tally_end = (struct contig_pos){ 0, 0 };
    tls.cur_pos = g_unset_pos;
    tls.ts = malloc(n_samples * sizeof(struct tally_stats));
    tls.refseqs = NULL;
    tls.n_refseqs = 0;
    unsigned s;
    for (s = 0; s != n_samples; ++s)
        tls.ts[s] = (struct tally_stats){
            .pbqt_hash = kh_init(pbqt_h),
            .indel_ct_hash = kh_init(indel_ct_h),
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

    tls.null_ct = (struct base_count){ { 0, 0, 0, 0 }, 0, 0 };

    unsigned pd = pseudo_depth;
    tls.refsam_ct[0] = (struct base_count){ { pd, 0, 0, 0 }, 0, pd };
    tls.refsam_ct[1] = (struct base_count){ { 0, pd, 0, 0 }, 0, pd };
    tls.refsam_ct[2] = (struct base_count){ { 0, 0, pd, 0 }, 0, pd };
    tls.refsam_ct[3] = (struct base_count){ { 0, 0, 0, pd }, 0, pd };
    tls.refsam_ct[4] = (struct base_count){ { 0, 0, 0, 0 }, 0, 0 };
}


void
batch_pileup_thread_free()
{
    fasta_thread_free();

    khiter_t it;
    for (it = kh_begin(tls.seq_hash); it != kh_end(tls.seq_hash); ++it)
        if (kh_exist(tls.seq_hash, it))
            free((char *)kh_key(tls.seq_hash, it));

    kh_destroy(str_h, tls.seq_hash);

    unsigned s;
    for (s = 0; s != tls.n_samples; ++s) {
        struct tally_stats *ts = &tls.ts[s];
        kh_destroy(pbqt_h, ts->pbqt_hash);

        khash_t(indel_ct_h) *ich = ts->indel_ct_hash;
        for (it = kh_begin(ich); it != kh_end(ich); ++it)
            if (kh_exist(ich, it)) {
                struct indel_count_ary *ia = &kh_val(ich, it);
                free(ia->i);
            }
        kh_destroy(indel_ct_h, ich);

        khash_t(olap_h) *oh = ts->overlap_hash;
        for (it = kh_begin(oh); it != kh_end(oh); ++it)
            if (kh_exist(oh, it)) {
                free((void *)kh_key(oh, it));
                free((void *)kh_val(oh, it));
            }
        kh_destroy(olap_h, oh);
        free(ts->base_ct);
        free(ts->bqs_ct);
        free(ts->indel_ct);
    }
}


static void
process_bam_block(char *rec, char *end,
                  struct contig_region *qbeg,
                  struct contig_region *qend,
                  struct contig_span loaded_span,
                  struct tally_stats *ts);


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

/* load specific ranges of reference sequence. [qbeg, qend) defines
   the total set of (non-overlapping) ranges to consider.  subset
   defines the overlapping intersection of these ranges that will be
   loaded into tls.refseqs.  assume that each interval in [qbeg, qend)
   is on one contig, but that 'subset' may span multiple contigs. */
void
pileup_load_refseq_ranges(struct bam_scanner_info *bsi)
{
    if (bsi->qbeg == bsi->qend) return;

    struct contig_region *q, *qlo, *qhi;
    find_intersecting_span(bsi->qbeg, bsi->qend, bsi->loaded_span, &qlo, &qhi);
    if (qlo == qhi) return;

    /* find the first region in [bsi->qbeg, bsi->qend) */
    unsigned r = 0, alloc = 0;
    for (q = qlo; q != qhi; ++q) {
        ALLOC_GROW(tls.refseqs, r + 1, alloc);
        tls.refseqs[r++] = (struct contig_fragment){ .reg = *q, .seq = NULL };
    }
    tls.n_refseqs = r;

    /* alter ranges of first and last (may be the same range) */
    struct contig_fragment *adj = &tls.refseqs[0];
    struct contig_pos new_beg = 
        MAX_CONTIG_POS(CONTIG_REGION_BEG(adj->reg), bsi->loaded_span.beg);
    assert(new_beg.tid == adj->reg.tid);
    adj->reg.beg = new_beg.pos;

    /* adjust last fragment */
    adj = &tls.refseqs[tls.n_refseqs - 1];

    struct contig_pos new_end = 
        MIN_CONTIG_POS(CONTIG_REGION_END(adj->reg), bsi->loaded_span.end);

    assert(new_end.tid == adj->reg.tid);
    adj->reg.end = new_end.pos;

    /* now load all sequences */
    for (r = 0; r != tls.n_refseqs; ++r) {
        struct contig_fragment *frag = &tls.refseqs[r];
        char *seq = fasta_fetch_iseq(frag->reg.tid,
                                     frag->reg.beg, 
                                     frag->reg.end);
        if ((frag->reg.end - frag->reg.beg) <= FRAGMENT_MAX_INLINE_SEQLEN) {
            strcpy(frag->buf, seq);
            frag->seq = frag->buf;
            free(seq);
        } else {
            frag->seq = seq;
        }
    }
    tls.cur_refseq = tls.refseqs;
}


/* release refseqs. */
static void
free_refseq_ranges()
{
    unsigned r;
    for (r = 0; r != tls.n_refseqs; ++r)
        if (tls.refseqs[r].seq != tls.refseqs[r].buf)
            free(tls.refseqs[r].seq);
    free(tls.refseqs);
    tls.refseqs = NULL;
    tls.n_refseqs = 0;
}


/* perform entire tally phase, for basecalls and for indels for one
   sample. update tls.tally_end to reflect furthest start position of
   any bam record seen. */
void
pileup_tally_stats(const struct managed_buf bam, 
                   struct bam_scanner_info *bsi,
                   unsigned s)
{
    process_bam_block(bam.buf, bam.buf + bam.size,
                      bsi->qbeg, bsi->qend, bsi->loaded_span,
                      &tls.ts[s]);
}    



    
/* return BAM bitflag-encoded form of current refbase. */
static unsigned
get_cur_refbase_code16()
{
    assert(tls.cur_pos.tid == tls.cur_refseq->reg.tid);
    unsigned offset = tls.cur_pos.pos - tls.cur_refseq->reg.beg;
    char refbase = tls.cur_refseq->seq[offset];
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
    if (s == REFERENCE_SAMPLE)
        return tls.refsam_ct[get_cur_refbase_code5()];

    struct tally_stats *ts = &tls.ts[s];
    assert(ts->base_cur == ts->base_end ||
           cmp_contig_pos(tls.cur_pos, ts->base_cur->cpos) < 1);

    if (ts->base_cur == ts->base_end
        || less_contig_pos(tls.cur_pos, ts->base_cur->cpos))
        return tls.null_ct;

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
        pk.k = kh_key(ph, it);
        pos = (struct contig_pos){ pk.v.tid, pk.v.pos };
        if (less_contig_pos(pos, tls.tally_end)) {
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
   pairs based on indel type. cts1 and cts2 are the indels from two
   samples at a given genomic position. */
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

    /* traverse both arrays in tandem, advancing one or both depending
       on comparison order.  */
    int cmp;
    while (ic0 != ie0 || ic1 != ie1) {
        cmp = (ic0 == ie0 ? -1 
               : (ic1 == ie1 ? 1 : less_indel(ic0->idl, ic1->idl)));
        switch (cmp) {
        case -1: 
            ip->count[0] = ic0->ct;
            ip->count[1] = 0; 
            ip->indel = ic0++->idl;
            break;
        case 1 : 
            ip->count[1] = ic1->ct;
            ip->count[0] = 0;
            ip->indel = ic1++->idl;
            break;
        case 0 :
            ip->count[0] = ic0->ct;
            ip->count[1] = ic1->ct; 
            ip->indel = ic0->idl; ++ic0; ++ic1;
            break;
        }
        ++ip;
    }            
    *n_pair_cts = ip - *pair_cts;
}


/* */
static inline char
bqs_count_to_call(struct bqs_count bc, unsigned refbase_i)
{
    static char match[] = ",.";
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


/* reify the indel, allocating and returning an indel_seq, using the
   current position.  Note: this function doesn't actually need the
   current position in order to retrieve sequences for insertions, but
   it does for deletions. */
struct indel_seq *
pileup_current_indel_seq(struct indel *idl)
{
    struct indel_seq *isq = malloc(sizeof(struct indel_seq) + idl->length + 1);
    isq->is_ins = idl->is_ins;
    if (isq->is_ins) {
        const char *ins_seq = kh_key(tls.seq_hash, idl->ins_itr);
        strcpy(isq->seq, ins_seq);
    } else {
        char *del_seq = fasta_fetch_iseq(tls.cur_pos.tid,
                                         tls.cur_pos.pos,
                                         tls.cur_pos.pos + idl->length);
        strcpy(isq->seq, del_seq);
        free(del_seq);
    }
    return isq;
}



/* produce pileup data (calls, quals, and read depths) from current
   position for sample s. manage buffer reallocation of pd fields. */
void
pileup_current_data(unsigned s, struct pileup_data *pd)
{
    if (s == REFERENCE_SAMPLE) {
        pd->calls.size = 3;
        ALLOC_GROW(pd->calls.buf, pd->calls.size, pd->calls.alloc);
        strncpy(pd->calls.buf, "REF", 3);

        pd->quals.size = 3;
        ALLOC_GROW(pd->quals.buf, pd->quals.size, pd->quals.alloc);
        strncpy(pd->quals.buf, "REF", 3);

        pd->n_match_lo_q = 0;
        pd->n_match_hi_q = pseudo_depth;
        pd->n_indel = 0;
        return;
    }        

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
    pd->n_indel = 0;
    for (b = 0; b != n_indel_ct; ++b) {
        indel_len_total += indel_ct[b].ct * (indel_ct[b].idl.length + 10);
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
    free(base_ct);

    /* print out all indel representations */
    struct indel_seq *isq;
    static char sign[] = "-+";
    for (b = 0; b != n_indel_ct; ++b) {
        isq = pileup_current_indel_seq(&indel_ct[b].idl);
        pd->calls.buf[p++] = sign[(int)isq->is_ins];
        p += sprintf(pd->calls.buf + p, "%Zu", strlen(isq->seq));
        p += strlen(isq->seq);
        free(isq);
    }
    pd->calls.size = p;
    free(indel_ct);
}


void
init_pileup_data(struct pileup_data *pd)
{
    pd->calls = (struct managed_buf){ NULL, 0, 0 };
    pd->quals = (struct managed_buf){ NULL, 0, 0 };
    pd->n_match_lo_q = 0;
    pd->n_match_hi_q = 0;
    pd->n_indel = 0;
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


static void
incr_indel_count_aux(struct tally_stats *ts,
                     struct contig_pos pos,
                     int len, 
                     uint8_t *bam_seq,
                     unsigned bam_rec_start);


static void
incr_insert_count(struct tally_stats *ts,
                  struct contig_pos pos,
                  unsigned len, 
                  uint8_t *bam_seq,
                  unsigned bam_rec_start)
{
    return incr_indel_count_aux(ts, pos, (int)len, bam_seq, bam_rec_start);
}


static void
incr_delete_count(struct tally_stats *ts,
                  struct contig_pos pos,
                  unsigned len)
{
    return incr_indel_count_aux(ts, pos, -(int)len, NULL, 0);
}


/* traverse b, tallying the match blocks into ts->pbqt_hash, and the
   indels into ts->indel_ct_hash and tls.seq_hash as necessary */
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
    
    /* initialize the constant parts of 'stat' here, then reuse
       this key in the loop */
    union pbqt_key stat = { .v = { .tid = tid, .strand = strand } };
    for (c = 0; c != b->core.n_cigar; ++c) {
        op = bam_cigar_op(cigar[c]);
        ln = bam_cigar_oplen(cigar[c]);

        khiter_t it;
        int32_t qe;
        if (op == BAM_CMATCH || op == BAM_CEQUAL || op == BAM_CDIFF) {
            qe = qpos + ln;
            for (q = qpos, r = rpos; q != qe; ++q, ++r) {
                stat.v.pos = r;
                stat.v.base = bam_seqi(bam_seq, q);
                stat.v.qual = bam_qual[q];
                if ((it = kh_get(pbqt_h, ph, stat.k)) == kh_end(ph)) {
                    it = kh_put(pbqt_h, ph, stat.k, &ret);
                    kh_val(ph, it) = 0;
                }
                kh_val(ph, it)++;
            }
        } else if (op == BAM_CINS) {
            struct contig_pos ins_pos = { tid, rpos };
            incr_insert_count(ts, ins_pos, ln, bam_seq, qpos);
        } else if (op == BAM_CDEL) {
            struct contig_pos del_pos = { tid, rpos };
            incr_delete_count(ts, del_pos, ln);
        } else
            ; /* All other operations (N, S, H, P) do not result in
                 tallying anything */
        
        /* but, we do advance the query */
        if (bam_cigar_type(op) & 1) qpos += ln;
        if (bam_cigar_type(op) & 2) rpos += ln;
    }
}




/* records are provided in [rec, end) ascending by start
   position. these are the minimal set of records that were identified
   by bam_scanner to overlap the region of interest defined by the
   intersection of 'loaded_span' and [qbeg, qend). but due to the
   limited resolution of the BSI index, there will be records in [rec,
   end) that do not overlap the region of interest. 

   loads each record, tallies its statistics if the record overlaps
   the region of interest. maintain a map of possibly overlapping read
   pairs and resolves their quality scores. return the upper bound
   position of complete tally statistics.  (this is the lowest start
   position of any records still in the overlaps_hash, or if empty,
   the start position of the last record parsed). */
static void
process_bam_block(char *rec, char *end,
                  struct contig_region *qbeg,
                  struct contig_region *qend,
                  struct contig_span loaded_span,
                  struct tally_stats *ts)
{
    /* the call to this function tallies all data in loaded_span.
       this signals the downstream consumers that it is complete to
       this point. */
    tls.tally_end = loaded_span.end;

    bam1_t b, b_mate;
    int ret;
    char *rec_next;

    struct contig_region *qcur, *qhi;
    find_intersecting_span(qbeg, qend, loaded_span, &qcur, &qhi);

    if (qcur == qhi) return;

    struct contig_pos rec_beg;
    khiter_t itr;
    while (rec != end) {
        rec_next = bam_parse(rec, &b);
        if (! rec_overlaps(&b, &qcur, qhi)) {
            if (qcur == qhi) break;
            else {
                rec = rec_next;
                continue;
            }
        }

        rec_beg = (struct contig_pos){ b.core.tid, b.core.pos };
        if (cmp_contig_pos(loaded_span.end, rec_beg) == -1)
            break; /* no more records overlapping the region of interest. done. */

        if (! (b.core.flag & BAM_FPAIRED) || ! (b.core.flag & BAM_FPROPER_PAIR)) {
            /* either this read is not paired or it is paired but not
               mapped in a proper pair.  tally and don't store. */
            process_bam_stats(&b, ts);
        } else { /* b is paired and mapped in proper pair */
            if (b.core.pos < b.core.mpos) { /* b is upstream mate */
                if (bam_endpos(&b) < b.core.mpos) { /* does not overlap with mate */
                    process_bam_stats(&b, ts);
                } else {
                    /* b overlaps mate.  since b is upstream, it
                       shouldn't be in overlaps hash. store; do not
                       tally yet. */
                    char *qname = strdup(bam_get_qname(&b));
                    itr = kh_put(olap_h, ts->overlap_hash, qname, &ret);

                    assert(ret == 1 || ret == 2);
                    kh_val(ts->overlap_hash, itr) = bam_duplicate_buf(rec);
                }
            } else {
                /* b is downstream mate. must search overlaps hash to
                   see if it overlaps with its upstream mate. */
                itr = kh_get(olap_h, ts->overlap_hash, bam_get_qname(&b));
                if (itr != kh_end(ts->overlap_hash)) {
                    (void)bam_parse(kh_val(ts->overlap_hash, itr), &b_mate);
                    // !!! NEED TO IMPLEMENT: tweak_overlap_quality(&b, &b_mate);
                    process_bam_stats(&b, ts);
                    process_bam_stats(&b_mate, ts);
                    free((void *)kh_key(ts->overlap_hash, itr));
                    free((void *)kh_val(ts->overlap_hash, itr));
                    kh_del(olap_h, ts->overlap_hash, itr);
                } else {
                    if (b.core.pos != b.core.mpos) {
                        /* downstream mate does not overlap with its
                           upstream partner.  or, upstream mate was
                           never loaded since it was not within the
                           region of interest. */
                        process_bam_stats(&b, ts);
                    } else {
                        /* b is first in the pair to be encountered.  store it */
                        char *qname = strdup(bam_get_qname(&b));
                        itr = kh_put(olap_h, ts->overlap_hash, qname, &ret);
                        assert(ret == 1 || ret == 2);
                        kh_val(ts->overlap_hash, itr) = bam_duplicate_buf(rec);
                    }
                }
            }
        }
        rec = rec_next;
    }
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
    
    ks_introsort(pbc_sort, i, ary);
    ts->base_ct = ary;
    ts->base_cur = ary;
    ts->base_end = ary + i;
}



static int
pos_iter_indel_count_less(struct pos_indel_count a, struct pos_indel_count b)
{
    return less_contig_pos(a.cpos, b.cpos)
        || (equal_contig_pos(a.cpos, b.cpos)
            && less_indel(a.ict.idl, b.ict.idl));
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
    unsigned alloc = kh_size(ih); /* lower limit estimate */
    ts->indel_ct = realloc(ts->indel_ct, alloc * sizeof(struct pos_indel_count));
    
    khiter_t it;
    unsigned i, a;
    union pos_key k;
    struct indel_count_ary ica;
    for (it = kh_begin(ih), i = 0; it != kh_end(ih); ++it)
        if (kh_exist(ih, it)) {
            k.k = kh_key(ih, it);
            if (less_contig_pos(k.v, tls.tally_end)) {
                ica = kh_val(ih, it);
                ALLOC_GROW(ts->indel_ct, i + ica.n, alloc);
                for (a = 0; a != ica.n; ++a)
                    ts->indel_ct[i++] = (struct pos_indel_count){ k.v, ica.i[a] };
            }
        }

    ks_introsort(pi_sort, i, ts->indel_ct);
    ts->indel_cur = ts->indel_ct;
    ts->indel_end = ts->indel_ct + i;
}


/* increment a count of a particular insertion in a sample s at
   pos. len is the length of the indel (negative means deletion,
   positive means insertion).  bam_seq is the seq field of the bam
   record provided by bam_get_seq(). if null, do not record the actual
   indel sequence.  */
static void
incr_indel_count_aux(struct tally_stats *ts,
                     struct contig_pos pos,
                     int len, 
                     uint8_t *bam_seq,
                     unsigned bam_rec_start)
{
    union pos_key k = { .v = pos };
    khash_t(indel_ct_h) *ih = ts->indel_ct_hash;

    khiter_t itr;
    int ret;
    if ((itr = kh_get(indel_ct_h, ih, k.k)) == kh_end(ih)) {
        itr = kh_put(indel_ct_h, ih, k.k, &ret);
        kh_val(ih, itr) = (struct indel_count_ary){ NULL, 0, 0 };
    }
    struct indel_count_ary ica = kh_val(ih, itr);

    char is_ins = (len >= 0);
    char *seq, *seqp;
    khiter_t seq_itr = 0; /* to suppress warnings */
    if (is_ins) {
        seq = malloc(len + 1);
        seqp = seq;
        unsigned q = bam_rec_start, qe = q + len;
        while (q != qe) {
            *seqp++ = seq_nt16_str[bam_seqi(bam_seq, q)]; /* macro!  don't touch q */
            ++q;
        }
        *seqp++ = '\0';
        if ((seq_itr = kh_get(str_h, tls.seq_hash, seq)) == kh_end(tls.seq_hash))
            seq_itr = kh_put(str_h, tls.seq_hash, seq, &ret);
        else
            free(seq);
    }

    /* linear search through the array */
    struct indel_count *ic = ica.i, *ice = ic + ica.n;
    unsigned ulen = abs(len);
    while (ic != ice) {
        if (ic->idl.is_ins != is_ins) continue;
        if (is_ins) {
            if (ic->idl.ins_itr == seq_itr) {
                ic->ct++;
                break;
            }
        } else {
            if (ic->idl.length == ulen) {
                ic->ct++;
                break;
            }
        }
    }            

    if (ic == ice) {
        /* didn't find one */
        ALLOC_GROW(ica.i, ica.n + 1, ica.m);
        if (is_ins)
            ica.i[ica.n] = (struct indel_count){ 
                { .is_ins = 1, .ins_itr = seq_itr },
                .ct = 1
            };
        else
            ica.i[ica.n] = (struct indel_count){
                { .is_ins = 0, .length = ulen },
                .ct = 1
            };
        ica.n++;
    }
}


void
pileup_clear_stats()
{
    unsigned s;
    khiter_t itr;
    for (s = 0; s != tls.n_samples; ++s) {
        struct tally_stats *ts = &tls.ts[s];
        kh_clear(pbqt_h, ts->pbqt_hash);

        /* clear finished entries in indel counts hash */
        khash_t(indel_ct_h) *ih = ts->indel_ct_hash;
        for (itr = kh_begin(ih); itr != kh_end(ih); ++itr) {
            if (! kh_exist(ih, itr)) continue;
            free(kh_val(ih, itr).i);
            kh_del(indel_ct_h, ih, itr);
        }
        kh_clear(indel_ct_h, ih);

        khash_t(olap_h) *oh = ts->overlap_hash;
        for (itr = kh_begin(oh); itr != kh_end(oh); ++itr)
            if (kh_exist(oh, itr)) {
                free((void *)kh_key(oh, itr));
                free((void *)kh_val(oh, itr));
                kh_del(olap_h, oh, itr);
            }

        /* completely clear temporary data */
        free(ts->base_ct);
        ts->base_ct = ts->base_cur = ts->base_end = NULL;
        free(ts->bqs_ct);
        ts->bqs_ct = ts->bqs_cur = ts->bqs_end = NULL;
        free(ts->indel_ct);
        ts->indel_ct = ts->indel_cur = ts->indel_end = NULL;
    }

    /* free refseqs */
    free_refseq_ranges();
    tls.cur_pos = g_unset_pos;
}



/* advance the internal position marker to the next available position
   in the region of interest. assume refseqs, cur_refseq and n_refseqs
   are initialized appropriately. */
static void
next_pos_aux()
{
    if (cmp_contig_pos(tls.cur_pos, g_unset_pos) == 0) {
        if (tls.n_refseqs)
            tls.cur_pos = CONTIG_REGION_BEG(tls.cur_refseq->reg);
        else
            tls.cur_pos = g_end_pos;
        return;
    }
    tls.cur_pos.pos++;
    if (cmp_contig_pos(tls.cur_pos, CONTIG_REGION_END(tls.cur_refseq->reg)) != -1) {
        ++tls.cur_refseq;
        if (tls.cur_refseq - tls.refseqs == tls.n_refseqs)
            tls.cur_pos = g_end_pos;
        else
            tls.cur_pos = CONTIG_REGION_BEG(tls.cur_refseq->reg);
    }
}


/* advance the internal position marker to the next valid position
   within the region of interest.  if the current position is already
   within it, do nothing.  for this purpose, consider the 'UNSET'
   position to be valid. return 0 if the position was okay, or 1 if it needed fixing. */
static int
fix_pos_aux()
{
    if (cmp_contig_pos(tls.cur_pos, g_end_pos) == 0)
        return 0;
    assert(tls.n_refseqs);
    if (cmp_contig_pos(tls.cur_pos, CONTIG_REGION_END(tls.cur_refseq->reg)) != -1) {
        next_pos_aux();
        return 1;
    }
    return 0;
}



/* update data iterators bqs_cur, base_cur, and indel_cur to point to
   the next available data at or beyond tls.cur_pos */
static void
update_data_iters()
{
    /* update all cur pointers. this works correctly even if  */
    unsigned s;
    struct tally_stats *ts;
    for (s = 0; s != tls.n_samples; ++s) {
        ts = &tls.ts[s];
        /* increment bqs iterator ('while' is used here, because bqs
           have multiple entries per position) */
        while (ts->bqs_cur != ts->bqs_end
               && less_contig_pos(ts->bqs_cur->cpos, tls.cur_pos))
            ++ts->bqs_cur;

        /* increment base iterator */
        while (ts->base_cur != ts->base_end
            && less_contig_pos(ts->base_cur->cpos, tls.cur_pos))
            ++ts->base_cur;

        /* increment indel iterator ('while' used here, because may be
           multiple entries per position) */
        while (ts->indel_cur != ts->indel_end
               && less_contig_pos(ts->indel_cur->cpos, tls.cur_pos))
            ++ts->indel_cur;
    }
}

/* advance internal position marker to the first position that has
   data (bqs, base or indel) in any sample.  if tls.cur_pos has data,
   this function does not update it. if there is no data altogether,
   tls.cur_pos is set to g_end_pos, indicating that we are at the end
   of data. this position may not be in the region of interest, and
   may need to be advanced further. assumes that bqs_cur, base_cur,
   and indel_cur are current w.r.t tls.cur_pos. */
static void
next_data_pos()
{
    struct contig_pos min_pos = g_end_pos;
    /* test bqs */
    unsigned s;
    struct tally_stats *ts;
    if (tls.ts[0].bqs_ct)
        for (s = 0; s != tls.n_samples; ++s) {
            ts = &tls.ts[s];
            if (ts->bqs_cur == ts->bqs_end) continue;
            min_pos = MIN_CONTIG_POS(min_pos, ts->bqs_cur->cpos);
            if (cmp_contig_pos(tls.cur_pos, min_pos) == 0)
                return;
        }

    /* test bases */
    if (tls.ts[0].base_ct)
        for (s = 0; s != tls.n_samples; ++s) {
            ts = &tls.ts[s];
            if (ts->base_cur == ts->base_end) continue;
            min_pos = MIN_CONTIG_POS(min_pos, ts->base_cur->cpos);
            if (cmp_contig_pos(tls.cur_pos, min_pos) == 0)
                return;
        }

    /* test indels */
    if (tls.ts[0].indel_ct)
        for (s = 0; s != tls.n_samples; ++s) {
            ts = &tls.ts[s];
            if (ts->indel_cur == ts->indel_end) continue;
            min_pos = MIN_CONTIG_POS(min_pos, ts->indel_cur->cpos);
            if (cmp_contig_pos(tls.cur_pos, min_pos) == 0)
                return;
        }

    /* the minimum position does not coincide with our current
       position, so update the current position. */
    tls.cur_pos = min_pos;
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
    next_pos_aux();

    if (skip_empty_loci) {
        unsigned was_fixed = 0;
        do {
            update_data_iters();
            next_data_pos();
            was_fixed = fix_pos_aux();
        } while (was_fixed);
    } else
        update_data_iters();

    return less_contig_pos(tls.cur_pos, g_end_pos);
}


void
pileup_current_info(struct pileup_locus_info *pli)
{
    const char *contig = fasta_get_contig(tls.cur_pos.tid);
    assert(contig != NULL);
    strcpy(pli->refname, contig);
    pli->refbase = seq_nt16_str[get_cur_refbase_code16()];
    pli->pos = tls.cur_pos.pos;
}


void
pileup_final_input()
{
    tls.tally_end = g_end_pos;
}
