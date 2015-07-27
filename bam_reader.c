#include "bam_reader.h"
#include "chunk_strategy.h"
#include "virtual_bound.h"

#include "sdg/kbtree.h"
#include "sdg/khash.h"
#include "htslib/sam.h"
#include "htslib/hts.h"
#include "htslib/hfile.h"
#include "defs.h"

#include <assert.h>

extern struct chunk_strategy cs_stats;

/* This is the conservative estimate, assuming the last block is the
   maximal size of 0x10000 */
static unsigned
tally_bgzf_bytes(hts_pair64_t *chunks, unsigned n_chunks)
{
    unsigned n, sz = 0;
    for (n = 0; n != n_chunks; ++n)
        sz += (chunks[n].v>>16) - (chunks[n].u>>16) + 0x10000;
    return sz;
}


/* Populates 'bins' array with the set of 'new bins' which is the
   difference between reg2bins(0, ti * 16384) and reg2bins(0, (ti-1) *
   16384).  Thus, ti is the zero-based index of a 16kb tiling
   window. Return the number of new bins found. */
static int
next_bins_aux(int ti, int *bins, int firstbin0, int n_small_bins)
{
    if (ti == n_small_bins) return 0;
    bins[0] = firstbin0 + ti;
    int i;
    for (i = 1; (ti % 8 == 1); ++i, ti >>= 3)
        bins[i] = hts_bin_parent(bins[i-1]);

    return i;
}

/* intervals are sorted in the usual way.  those which have a full
   containment relationship are deemed equal.*/
static int
cmp_ival(hts_pair64_t a, hts_pair64_t b)
{
    if (a.v <= b.u) return -1;
    if (b.v <= a.u) return 1;
    if (a.u < b.u && a.v < b.v) return -1;
    if (b.u < a.u && b.v < a.v) return 1;
    return 0;
}

KBTREE_INIT(itree, hts_pair64_t, cmp_ival);

/* add an interval to the itree, resolving overlaps. return the total
   interval length gained. */
static uint64_t
add_chunk_to_itree(kbtree_t(itree) *itree, hts_pair64_t q)
{
    hts_pair64_t *lw, *hi, *pk[2];
    uint64_t ins_sz = 0, del_sz = 0;
    int i, cont = 0;
    while (1) {
        kb_interval(itree, itree, q, &lw, &hi);
        pk[0] = lw; pk[1] = lw != hi ? hi : NULL;

        /* delete lw or hi if either are contained by q */
        cont = 0;
        for (i = 0; i != 2; ++i)
            if (pk[i] && q.u <= pk[i]->u && pk[i]->v <= q.v) {
                del_sz += pk[i]->v - pk[i]->u;
                kb_delp(itree, itree, pk[i]);
                cont = 1;
            }
        if (cont) continue;

        /* q doesn't contain lw or hi. edit q to eliminate overlaps */
        if (lw && q.u < lw->v) q.u = lw->v;
        if (hi && hi->u < q.v) q.v = hi->u;
        ins_sz = q.v - q.u;

        /* q does not overlap any neighbors */
        if (q.u < q.v) kb_put(itree, itree, q);
        break;
    }
    assert(ins_sz >= del_sz);
    return ins_sz - del_sz;
}

/* taken from htslib/hts.c */
typedef struct {
    int32_t m, n;
    uint64_t loff;
    hts_pair64_t *list;
} bins_t;

KHASH_MAP_INIT_INT(bin, bins_t)
typedef khash_t(bin) bidx_t;

typedef struct {
    int32_t n, m;
    uint64_t *offset;
} lidx_t;

struct __hts_idx_t {
    int fmt, min_shift, n_lvls, n_bins;
    uint32_t l_meta;
    int32_t n, m;
    uint64_t n_no_coor;
    bidx_t **bidx;
    lidx_t *lidx;
    uint8_t *meta;
    struct {
        uint32_t last_bin, save_bin;
        int last_coor, last_tid, save_tid, finished;
        uint64_t last_off, save_off;
        uint64_t off_beg, off_end;
        uint64_t n_mapped, n_unmapped;
    } z; // keep internal states
};
/* end of section from htslib/hts.c */

/* taken from htslib/bgzf.c */
#define BLOCK_HEADER_LENGTH 18


/* add the set of chunk intervals from bin into the itree. resolve
   overlaps, tally the total length of intervals in the tree. */
static unsigned
add_bin_chunks_aux(int bin,
                   bidx_t *bidx, 
                   uint64_t min_off, 
                   kbtree_t(itree) *itree)
{
    khiter_t k;
    unsigned j, bgzf_bytes = 0;
    if ((k = kh_get(bin, bidx, bin)) != kh_end(bidx)) {
        bins_t *p = &kh_val(bidx, k);
        for (j = 0; j != p->n; ++j)
            if (p->list[j].v > min_off)
                bgzf_bytes += add_chunk_to_itree(itree, p->list[j]);
    }
    return bgzf_bytes;
}

/* same as reg2bins, except that the [tbeg, tend) subset of contiguous
   16kb windows is given.  windows 0, 1, and 2 are respectively [0,
   16kb), [16kb, 32kb), [32kb, 48kb).  n_lvls: number of levels.
   i.e. 6 levels gives levels [0, 6) */
static inline int
win2bins(int tbeg, int tend, hts_itr_t *itr, int n_lvls)
{
    int l, t, max_lvl = n_lvls - 1, s = (max_lvl<<1) + max_lvl;
    if (tbeg == tend) return 0;
    assert(tbeg <= tend);
    assert(tend <= 1<<s);
    
    
    for (l = 0, t = 0; l != n_lvls; s -= 3, t += 1<<((l<<1)+l), ++l) {
        int b, e, n, i;
        b = t + (tbeg>>s); e = t + (tend>>s); n = e - b + 1;
        if (itr->bins.n + n > itr->bins.m) {
            itr->bins.m = itr->bins.n + n;
            kroundup32(itr->bins.m);
            itr->bins.a = (int*)realloc(itr->bins.a, sizeof(int) * itr->bins.m);
        }
        for (i = b; i <= e; ++i) itr->bins.a[itr->bins.n++] = i;
    }
    return itr->bins.n;
    
}


/* compute min_off. ripped from htslib/hts.c.  traverse bins in
   order of previous siblings, then parents (which is the order of
   descending loff) until a non-empty bin is found and loff is
   initialized. */
static uint64_t
find_min_offset(int tid, int ti, hts_idx_t *idx)
{
    bidx_t *bidx = idx->bidx[tid];
    int bin = hts_bin_first(idx->n_lvls) + ti;
    khiter_t k;
    do {
        int first;
        k = kh_get(bin, bidx, bin);
        if (k != kh_end(bidx)) break;
        first = (hts_bin_parent(bin)<<3) + 1;
        if (bin > first) --bin;
        else bin = hts_bin_parent(bin);
    } while (bin);
    if (bin == 0) k = kh_get(bin, bidx, bin);
    return k != kh_end(bidx) ? kh_val(bidx, k).loff : 0;
}


/* compute initial set of bins, add those chunks to the tree. return
   the total number of bytes added. */
static unsigned
compute_initial_bins(int tid, 
                     int ti,
                     uint64_t min_off, 
                     kbtree_t(itree) *itree, 
                     hts_idx_t *idx)
{
    hts_itr_t *itr = calloc(1, sizeof(hts_itr_t));

    int actual_n_lvls = idx->n_lvls + 1; /* semantic correction */
    win2bins(ti, ti + 1, itr, actual_n_lvls);
    unsigned bgzf_bytes = 0;

    /* add initial chunks for [ti, ti + 1) */
    int i;
    for (i = 0; i != itr->bins.n; ++i)
        bgzf_bytes += add_bin_chunks_aux(itr->bins.a[i], 
                                         idx->bidx[tid], min_off, itree);
    free(itr);
    return bgzf_bytes;
}


/* transfer all chunks in itree to chunks, appending from initial
   position *n_chunks and updating both chunks and n_chunks. */
static void
accumulate_chunks(kbtree_t(itree) *itree,
                  hts_pair64_t **chunks,
                  unsigned *n_chunks)
{
    unsigned beg = *n_chunks;
    *n_chunks += kb_size(itree);
    *chunks = realloc(*chunks, *n_chunks * sizeof((*chunks)[0]));

    kbitr_t tree_itr;
    kb_itr_first_itree(itree, &tree_itr);
    int iret = 1; /* */
    unsigned i;
    for (i = beg; i != *n_chunks; ++i) {
        assert(iret != 0);
        (*chunks)[i] = kb_itr_key(hts_pair64_t, &tree_itr);
        iret = kb_itr_next_itree(itree, &tree_itr);
    }
}

struct virt_less_range_par {
    struct pair_ordering_range *ary;
    struct pair_ordering q;
};

/* uses par both as a source of elements elem and a query element q.
   return 1 if q < elem[pos], 0 otherwise.
 */
static int
less_rng_end(unsigned pos, void *par)
{
    struct virt_less_range_par *vl = par;
    return cmp_pair_ordering(&vl->q, &vl->ary[pos].end) < 0;
}

#define MIN_PAIR_ORD(a, b) (cmp_pair_ordering(&(a), &(b)) < 0 ? (a) : (b))
#define MAX_PAIR_ORD(a, b) (cmp_pair_ordering(&(a), &(b)) < 0 ? (b) : (a))

#define MIN(a, b) ((a) < (b) ? (a) : (b))

/* using forward scanning through the index, return the largest
   position end such that the total size of bgzf blocks overlapping
   [beg, end) is at or slightly above min_wanted_bytes. use [qbeg,
   qend) to define the subset of ranges to consider.  max_end defines
   the maximum logical end range to consider.

   (*chunks)[c] = chunk c
   (*n_chunks) = number of chunks
*/
static struct pair_ordering
hts_size_to_range(struct pair_ordering beg,
                  struct pair_ordering max_end,
                  unsigned min_wanted_bytes,
                  struct pair_ordering_range *qbeg,
                  struct pair_ordering_range *qend,
                  hts_idx_t *idx,
                  hts_pair64_t **chunks,
                  unsigned *n_chunks)
{
    /* constants */
    int firstbin0 = hts_bin_first(idx->n_lvls);
    int firstbin1 = hts_bin_first(idx->n_lvls - 1);
    int n_small_bins = firstbin0 - firstbin1;
    int n_bins, *bins = malloc(idx->n_lvls * sizeof(int));

    /* initialize qcur */
    struct virt_less_range_par vpar = { qbeg, beg };
    struct pair_ordering_range *qcur = qbeg
        + virtual_upper_bound(0, qend - qbeg, less_rng_end, &vpar);
        
    int ms = idx->min_shift;
    struct pair_ordering cpos;
    struct pair_ordering_range crng; /* current range */
    unsigned bgzf_bytes;

    for (*n_chunks = 0, bgzf_bytes = 0; qcur != qend; ++qcur) {
        crng = (struct pair_ordering_range){ MAX_PAIR_ORD(qcur->beg, beg), qcur->end };
        int ti, b, tid, tid_end = qcur->end.hi + 1;
        for (tid = qcur->beg.hi; tid != tid_end; ++tid) {
            int ti_beg = tid == crng.beg.hi ? crng.beg.lo>>ms : 0;
            int ti_end = tid == crng.end.hi ? crng.end.lo>>ms : n_small_bins;
            if (ti_beg == ti_end) continue;

            uint64_t min_off = find_min_offset(tid, ti_beg, idx);

            kbtree_t(itree) *itree = kb_init(itree, 128);
            bgzf_bytes += compute_initial_bins(tid, ti_beg, min_off, itree, idx);

            cpos.hi = tid;
            cpos.lo = (ti_beg + 1)<<ms;

            for (ti = ti_beg + 1; ti <= ti_end; ++ti) {
                if (bgzf_bytes >= min_wanted_bytes
                    || cmp_pair_ordering(&cpos, &max_end) >= 0) {
                    /* break if window touches or exceeds max_end, or
                       we have enough content */
                    accumulate_chunks(itree, chunks, n_chunks);
                    kb_destroy(itree, itree);
                    goto END;
                }

                /* process this 16kb window */
                n_bins = next_bins_aux(ti, bins, firstbin0, n_small_bins);
                for (b = 0; b != n_bins; ++b)
                    bgzf_bytes += add_bin_chunks_aux(bins[b], idx->bidx[tid], min_off, itree);
                cpos.lo = (ti + 1)<<ms;
            }
            accumulate_chunks(itree, chunks, n_chunks);
            kb_destroy(itree, itree);
        }
        /* beg is only relevant within the first value of qcur */
        beg = (struct pair_ordering){ 0, 0 }; 
    }
 END:
    free(bins);
    return cpos;
}

/* taken from htslib/bgzf.c */
static inline int
unpackInt16(const uint8_t *buffer)
{
    return buffer[0] | buffer[1] << 8;
}


/* scans the bam index to find all chunks covering a range starting at
   cs_stats.pos.  scans forward until *close to* but less than
   max_bytes are found. stores the found chunks in par's fields, and
   updates cs_stats.pos to the newest position.  */
#define STEP_DOWN 0.975
void
bam_scanner(void *par, unsigned max_bytes)
{
    struct bam_reader_par *bp = par;
    float mul = STEP_DOWN;
    unsigned tmp_bytes, min_actual_bytes = UINT_MAX, min_wanted_bytes;
    struct pair_ordering tmp_end, min_end = { SIZE_MAX, SIZE_MAX };
    do {
        min_wanted_bytes = max_bytes * mul;
        unsigned s;
        for (s = 0; s != bp->n; ++s) {
            tmp_end =
                hts_size_to_range(cs_stats.pos,
                                  min_end,
                                  min_wanted_bytes,
                                  bp->qbeg,
                                  bp->qend,
                                  bp->m[s].idx, 
                                  &bp->m[s].chunks,
                                  &bp->m[s].n_chunks);
            min_end = MIN_PAIR_ORD(min_end, tmp_end);
            tmp_bytes = tally_bgzf_bytes(bp->m[s].chunks, bp->m[s].n_chunks);
            min_actual_bytes = MIN(min_actual_bytes, tmp_bytes);
        }
        mul *= STEP_DOWN;
    } while (min_actual_bytes > max_bytes);

    cs_stats.pos = min_end;
}


/* reads the virtual offset ranges serialized in buf.  returns a new
   buffer holding the voffset pairs.  
   bytes.  allocates new memory and points *chunks to it. */
static hts_pair64_t *
read_voffset_chunks(const struct managed_buf *buf,
                    unsigned *n_chunks,
                    size_t *n_bytes_consumed)
{
    *n_chunks = *(unsigned *)buf->buf;
    unsigned
        sz1 = sizeof(*n_chunks), 
        sz2 = *n_chunks * sizeof(hts_pair64_t);
    hts_pair64_t *chunks = malloc(sz2);
    memcpy(chunks, buf->buf + sz1, sz2);
    *n_bytes_consumed = sz1 + sz2;
    return chunks;
}


/* serialize chunks into buf */
static void
write_voffset_chunks(hts_pair64_t *chunks,
                     unsigned n_chunks,
                     struct managed_buf *buf)
{
    buf->size = sizeof(unsigned) + n_chunks * sizeof(hts_pair64_t);
    ALLOC_GROW(buf->buf, buf->size, buf->alloc);
    *(unsigned *)buf->buf = n_chunks;
    size_t sz = n_chunks * sizeof(chunks[0]);
    memcpy(buf->buf + sizeof(unsigned), chunks, sz);
}


/* bgzf points to the start of a bgzf block.  reads the block size
   field of this record, per the BAM spec. */
static inline int
bgzf_block_size(char *bgzf)
{
    uint8_t *buf = (uint8_t *)bgzf + 16;
    return unpackInt16(buf) + 1;
}


/* read recorded compressed blocks into bufs. plus space for the
   initial voffset information.  Warning: Does not check that hread
   calls return the correct number of characters read.  */
void
bam_reader(void *par, struct managed_buf *bufs)
{
    struct bam_reader_par *bp = par;
    unsigned s, c;
    int block_span, bgzf_block_len, remain;
    ssize_t n_chars_read;
    char *wp;

    unsigned bytes_bound, chunk_hdr_sz;

    for (s = 0; s != bp->n; ++s) {
        bytes_bound = tally_bgzf_bytes(bp->m[s].chunks, bp->m[s].n_chunks);
        chunk_hdr_sz = sizeof(unsigned) + bp->m[s].n_chunks * sizeof(bp->m[s].chunks[0]);

        ALLOC_GROW(bufs[s].buf, bytes_bound + chunk_hdr_sz, bufs[s].alloc);
        
        write_voffset_chunks(bp->m[s].chunks, bp->m[s].n_chunks, &bufs[s]);
        wp = bufs[s].buf + bufs[s].size;

        hFILE *fp = bp->m[s].bgzf->fp;

        for (c = 0; c != bp->m[s].n_chunks; ++c) {
            bgzf_seek(bp->m[s].bgzf, bp->m[s].chunks[c].u, SEEK_SET);

            /* read all but last block in this chunk */
            block_span = (bp->m[s].chunks[c].v>>16) - (bp->m[s].chunks[c].u>>16);
            n_chars_read = hread(fp, wp, block_span);
            
            /* determine size of last block and read remainder */
            wp += block_span;
            n_chars_read = hread(fp, wp, BLOCK_HEADER_LENGTH);

            bgzf_block_len = bgzf_block_size(wp);
            wp += BLOCK_HEADER_LENGTH;
            remain = bgzf_block_len - BLOCK_HEADER_LENGTH;

            n_chars_read = hread(fp, wp, remain);
            assert(n_chars_read == remain);

            wp += remain;
        }
        bufs[s].size = wp - bufs[s].buf;
    }
}

static size_t
inflate_bgzf_block(char *in, int block_length, char *out);


/* bgzf contains first a header consisting of <n_blocks> and then a
   number of hts_pair64_t virtual offset pairs.  following that is
   actual bgzf blocks to be inflated.  stores inflated BAM records in
   bam.  manage the size of bam.  only copy the portions of blocks
   defined by the virtual offsets. */
void
bam_inflate(const struct managed_buf *bgzf,
            struct managed_buf *bam)
{
    unsigned c;
    int64_t c_beg, c_end; /* sub-regions in a chunk that we want */
    size_t n_copy;
    uint64_t bs; /* compressed block size */

    uint64_t coff, coff_beg, coff_end; /* current and end coffset as we traverse
                                          the chunks within our block ranges */
    
    /* parse the first part of the bgzf buffer to get the set of
       chunks */
    unsigned n_chunks;
    size_t n_bytes_consumed;
    hts_pair64_t *chunks = 
        read_voffset_chunks(bgzf, &n_chunks, &n_bytes_consumed);

    bam->size = 0;
    char *in = bgzf->buf + n_bytes_consumed, *out = bam->buf;
    size_t n_bytes;

    for (c = 0; c != n_chunks; ++c) {
        coff_beg = chunks[c].u>>16;
        coff_end = chunks[c].v>>16; /* offset of the beginning of the last block in this chunk */
        coff = coff_beg;
        do {
            ALLOC_GROW_REMAP(bam->buf, out, bam->size + BGZF_MAX_BLOCK_SIZE, bam->alloc);
            bs = bgzf_block_size(in);

            n_bytes = inflate_bgzf_block(in, bs, out);
            c_beg = (coff == coff_beg ? chunks[c].u & 0xffff : 0);
            c_end = (coff == coff_end ? chunks[c].v & 0xffff : n_bytes);
            n_copy = c_end - c_beg; /* number of bytes to copy to output buffer */

            if (c_beg != 0) memmove(out, out + c_beg, n_copy);

            bam->size += n_copy;
            in += bs;
            out += n_copy;
            coff += bs;
        } while (coff <= coff_end);
    }
}


static size_t
inflate_bgzf_block(char *in, int block_length, char *out)
{
    z_stream zs;
    zs.zalloc = NULL;
    zs.zfree = NULL;
    zs.next_in = (Bytef*)in + 18;
    zs.avail_in = block_length - 16;
    zs.next_out = (Bytef*)out;
    zs.avail_out = BGZF_MAX_BLOCK_SIZE;

    if (inflateInit2(&zs, -15) != Z_OK)
        return -1;

    if (inflate(&zs, Z_FINISH) != Z_STREAM_END) {
        inflateEnd(&zs);
        return -1;
    }
    if (inflateEnd(&zs) != Z_OK) {
        return -1;
    }
    return zs.total_out;
}


/* parse a raw record into b. return the position of the next record
   in the block */
char *
bam_parse(char *raw, bam1_t *b)
{
    bam1_core_t *c = &b->core;
    b->data = (uint8_t *)raw + 36; /* do not allocate data, just point to raw buffer. */
    int32_t block_len = *(int32_t *)raw;
    

#ifdef DEP_BIG_ENDIAN
    uint32_t x[8];
    ed_swap_4p(&block_len);
    unsigned i;
    memcpy(x, raw + 4, 32);
    for (i = 0; i != 8; ++i) ed_swap_4p(x + i);
    swap_data(&c, block_len - 32, raw, 0);
#else
    uint32_t *x = raw + 4;
#endif

    c->tid = x[0]; c->pos = x[1];
    c->bin = x[2]>>16; c->qual = x[2]>>8&0xff; c->l_qname = x[2]&0xff;
    c->flag = x[3]>>16; c->n_cigar = x[3]&0xffff;
    c->l_qseq = x[4];
    c->mtid = x[5]; c->mpos = x[6]; c->isize = x[7];

    b->l_data = b->m_data = block_len - 32;
    return raw + block_len + 4;
}


/* initialize bs from bam_file.  also, expects <bam_file>.bai to exist
   as well */
void
bam_stats_init(const char *bam_file, struct bam_stats *bs)
{
    bs->bgzf = bgzf_open(bam_file, "r");
    char *bam_index_file = malloc(strlen(bam_file) + 5);
    strcpy(bam_index_file, bam_file);
    strcat(bam_index_file, ".bai");

    bs->idx = bam_index_load(bam_index_file);
    bs->hdr = bam_hdr_read(bs->bgzf);
    bs->chunks = NULL;
    bs->n_chunks = 0;
}


void
bam_stats_free(struct bam_stats *bs)
{
    hts_idx_destroy(bs->idx);
    bam_hdr_destroy(bs->hdr);
    bgzf_close(bs->bgzf);
    free(bs->chunks);
}


#if 0
void bam_reader(void *par, struct managed_buf *bufs)
{
    struct bam_reader_par *bp = par;
    unsigned s, i;
    char *bgzf_buf = NULL, *bgzf_end, *rec_ptr;
    uint64_t n_bgzf_bytes;
    uint8_t header[BLOCK_HEADER_LENGTH];
    int bgz_block_len, rec_block_len;
    for (s = 0; s != bp->n_idx; ++s) {
        /* point fp->compressed_block field to populated raw
           buffer. */
        cb_old = fp->compressed_block;
        rec_ptr = bufs[s].buf;
        for (i = 0; i != it.n_off; ++i) {
            block_span = (it.off[i].v>>16) - (it.off[i].u>>16);
            block_offset = it.off[i].u & 0xFFFF;
            is_last_block = 0;
            n_bgzf_bytes =  block_span + 0x10000;
            bgzf_buf = realloc(bgzf_buf, n_bgzf_bytes);
            bgzf_seek(fp, it.off[i].u, SEEK_SET);
            hread(fp->fp, bgzf_buf, n_bgzf_bytes);
            fp->compressed_block = bgzf_buf;
            
            while (! is_last_block) {
                bgzf_block_len = unpackInt16((uint8_t *)&fp->compressed_block[16]) + 1;
                if (rec_block_len = inflate_block(fp->fp, bgzf_block_len) < 0) return -1;
                fp->compressed_block += bgzf_block_len;
                is_last_block = (fp->compressed_block - bgzf_buf) == block_span;
                block_offset = is_last_block ? it.off[i].v & 0xFFFF : 0;
                
                memcpy(rec_ptr, 
                       fp->uncompressed_block + block_offset,
                       rec_block_len - block_offset);
                
                rec_ptr += rec_block_len - block_offset;
            }
        }
        fp->compressed_block = cb_old;
        /* is fp in a consistent state now? */
    }

    free(bgzf_buf);
}
#endif
