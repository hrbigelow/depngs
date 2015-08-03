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

/* extract compressed offset from voffset */
#define BAM_COFFSET(v) ((v)>>16)

/* extract uncompressed offset from voffset */
#define BAM_UOFFSET(v) ((v) & 0xFFFF)

#define MIN(a, b) ((a) < (b) ? (a) : (b))
#define MAX(a, b) ((a) < (b) ? (b) : (a))

#define MIN_CONTIG_POS(a, b) (cmp_contig_pos((a), (b)) < 0 ? (a) : (b))
#define MAX_CONTIG_POS(a, b) (cmp_contig_pos((a), (b)) < 0 ? (b) : (a))

/* Estimate the total size in bytes of blocks in this chunk. If the
   within-block offset of the last block is zero, then we don't need
   to parse that last block.  Otherwise, we need to parse it, and
   cannot know its size in advance, so use 0x10000 as the maximum
   possible size. */
static unsigned
chunk_bytes(hts_pair64_t chunk)
{
    return 
        BAM_COFFSET(chunk.v) 
        - BAM_COFFSET(chunk.u) 
        + (BAM_UOFFSET(chunk.v) ? BGZF_MAX_BLOCK_SIZE : 0);
}

/* Return an estimate for the bytes required to hold the BGZF blocks
   contained in the chunks.  Each chunk
   assuming the last block is the maximal size of 0x10000 */
static unsigned
tally_bgzf_bytes(hts_pair64_t *chunks, unsigned n_chunks)
{
    unsigned c, sz = 0;
    for (c = 0; c != n_chunks; ++c)
        sz += chunk_bytes(chunks[c]);
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

/* Return 0 if chunks a and b share one or more BGZF blocks.  Return
   -1 if chunk a's blocks are all earlier in the file than chunk b's
   blocks, and 1 if the reverse is true.

   Note.  It may happen that a.v < b.u, but BAM_COFFSET(a.v) ==
   BAM_COFFSET(b.u).  If this is the case, and if BAM_UOFFSET(a.v) != 0, then
   this means that chunk a and chunk b both contain the block at
   BAM_COFFSET(a.v).  However, if BAM_UOFFSET(a.v) == 0, then chunk a does NOT
   contain the block at BAM_COFFSET(a.v) and so chunks a and b are deemed
   not to overlap.  */
static int
cmp_bgzf_chunks(hts_pair64_t a, hts_pair64_t b)
{
    if (BAM_COFFSET(a.v) < BAM_COFFSET(b.u) ||
        (BAM_COFFSET(a.v) == BAM_COFFSET(b.u) && BAM_UOFFSET(a.v) == 0))
        return -1;
    if (BAM_COFFSET(b.v) < BAM_COFFSET(a.u) ||
        (BAM_COFFSET(b.v) == BAM_COFFSET(a.u) && BAM_UOFFSET(b.v) == 0))
        return 1;

    return 0;
}


/* return 1 if c_query is 'contained' by c_target: that is, if the set
   of c_query's BGZF blocks is a subset of c_target's set of BGZF
   blocks. */
static int
bgzf_contains(hts_pair64_t c_target, hts_pair64_t c_query)
{
    return 
        BAM_COFFSET(c_target.u) <= BAM_COFFSET(c_query.u)
        && (BAM_COFFSET(c_query.v) < BAM_COFFSET(c_target.v)
            || (BAM_COFFSET(c_query.v) == BAM_COFFSET(c_target.v)
                && BAM_UOFFSET(c_target.v) != 0));
}


/* produce a new chunk (voffset pair) that equals the union of BGZF
   blocks defined by chunks c1 and c2.  The starting and ending
   uoffsets will be preserved.  However, total set of BAM records in
   the merged chunk may be larger, in case that the two chunks shared
   a block, but used disjoint subsets of records within that block. 

   Assumes that c1 and c2 intersect.  */
static hts_pair64_t
merge_chunks(hts_pair64_t c1, hts_pair64_t c2)
{
    return (hts_pair64_t){ .u = MIN(c1.u, c2.u), .v = MAX(c1.v, c2.v) };
}


KBTREE_INIT(itree_t, hts_pair64_t, cmp_bgzf_chunks);

/* add an interval to the itree, resolving overlaps. return the total
   size of blocks.

   Approach: query the itree for any item that compares equal to
   q. For any existing item e found, determine whether e contains q.
   If so, do nothing.  Otherwise, remove e, update del_sz, set q to
   merged(q, e) and repeat.  This procedure maintains the invariant
   that no two items in the tree compare equal (according to
   cmp_bgzf_chunks, which deems overlap to be equality).  Return the
   change in total size of all chunks in the itree (as measured by
   chunk_bytes()).  Note that the change may be negative!  This is
   because items may be merged by the addition of a new item, and
   chunk_bytes() provides a conservative over-estimate of size
   (assuming the last BGZF block is BGZF_MAX_BLOCK_SIZE) */
static int64_t
add_chunk_to_itree(kbtree_t(itree_t) *itree, hts_pair64_t q)
{
    int64_t ins_sz = 0, del_sz = 0;
    hts_pair64_t *ex; /* existing element */
    while (1) {
        ex = kb_get(itree_t, itree, q);
        if (ex) {
            if (bgzf_contains(*ex, q))
                break;
            else {
                del_sz += chunk_bytes(*ex);
                q = merge_chunks(*ex, q);
                kb_del(itree_t, itree, *ex);
                continue;
            }
        } else {
            ins_sz += chunk_bytes(q);
            kb_put(itree_t, itree, q);
            break;
        }
    }
    return ins_sz - del_sz;
}


static unsigned
total_itree_bytes(kbtree_t(itree_t) *itree)
{
    unsigned sz = 0;
    kbitr_t itr;
    int iret = 1;
    hts_pair64_t ck;
    kb_itr_first_itree_t(itree, &itr);
    while (iret) {
        ck = kb_itr_key(hts_pair64_t, &itr);
        sz += chunk_bytes(ck);
        iret = kb_itr_next_itree_t(itree, &itr);
    }
    return sz;
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
                   kbtree_t(itree_t) *itree)
{
    khiter_t k;
    unsigned j, bgzf_bytes = 0;
    if ((k = kh_get(bin, bidx, bin)) != kh_end(bidx)) {
        bins_t *p = &kh_val(bidx, k);
        for (j = 0; j != p->n; ++j)
            if (p->list[j].v > min_off) {
                int64_t before_bytes = total_itree_bytes(itree);
                int64_t bytes_added = add_chunk_to_itree(itree, p->list[j]);
                bgzf_bytes += bytes_added;
                int64_t after_bytes = total_itree_bytes(itree);
                assert(bytes_added == (after_bytes - before_bytes));
            }
    }
    return bgzf_bytes;
}


/* return v rounded to nearest multiple of 8^p. Hacker's Delight, 2nd
   ed., page 59 */
static inline unsigned
round_up_pow8(unsigned v, unsigned p)
{
    assert(p < 6);
    /* r1[k] = 1<<(3*k) - 1*/
    static unsigned r1[] = { 0, 7, 63, 511, 4095, 32767 };
    static unsigned r2[] = { ~0, ~7, ~63, ~511, ~4095, ~32767 };
    return (v + r1[p]) & r2[p];
}


/* returns an allocated array of the set of bins that overlap the
   contiguous tiling window range [tbeg, tend).  sets *n_bins to the
   number of bins found.  same as reg2bins, except that the [tbeg,
   tend) subset of contiguous 16kb windows is given.  windows 0, 1,
   and 2 are respectively [0, 16kb), [16kb, 32kb), [32kb, 48kb).
   n_lvls: number of levels.  i.e. 6 levels gives levels [0, 6).  */
static inline int *
win2bins(int tbeg, int tend, int n_lvls, unsigned *n_bins)
{
    if (tbeg == tend) {
        *n_bins = 0;
        return NULL;
    }

    int l, t, max_lvl = n_lvls - 1, s = (max_lvl<<1) + max_lvl;
    assert(tbeg <= tend);
    assert(tend <= 1<<s);

    unsigned n = 0, alloc = n_lvls;
    int *bins = malloc(alloc * sizeof(int));

    /* identify the bin range [b, e) at each level l such that each
       bin in the range overlaps at least one tiling window in the
       range [tbeg, tend). t gives the bin number of the first bin on
       level l. */
    for (l = 0, t = 0; l != n_lvls; s -= 3, t += 1<<((l<<1)+l), ++l) {
        int b, e, add, i;
        unsigned p = n_lvls - l - 1;
        b = t + (round_up_pow8(tbeg, p)>>s);
        e = t + (round_up_pow8(tend, p)>>s);
        add = e - b;
        ALLOC_GROW(bins, n + add, alloc);
        for (i = b; i != e; ++i) bins[n++] = i;
    }
    *n_bins = n;
    return bins;
}


/* compute min_off. tid is the zero-based index of the contig.  ti is
   the zero-based index of the 16kb tiling window.  finds the minimum
   voffset needed in order to search the BAM file for records
   overlapping this 16kb tiling window.  is ripped from htslib/hts.c.
   traverse bins in order of previous siblings, then parents (which is
   the order of descending loff) until a non-empty bin is found and
   loff is initialized. */
static uint64_t
find_min_offset(int tid, int ti, hts_idx_t *idx)
{
    bidx_t *bidx = idx->bidx[tid];
    int bin = hts_bin_first(idx->n_lvls) + ti;
    khiter_t itr;
    do {
        int first;
        itr = kh_get(bin, bidx, bin);
        if (itr != kh_end(bidx)) break;
        first = (hts_bin_parent(bin)<<3) + 1;
        if (bin > first) --bin;
        else bin = hts_bin_parent(bin);
    } while (bin);
    if (bin == 0) itr = kh_get(bin, bidx, bin);
    return itr != kh_end(bidx) ? kh_val(bidx, itr).loff : 0;
}


/* compute initial set of bins that have content in tiling windows
   [ti, ti + 1) (?), add those chunks to the tree. return the total
   number of bytes added. */
static unsigned
compute_initial_bins(int tid, 
                     int ti,
                     uint64_t min_off, 
                     kbtree_t(itree_t) *itree, 
                     hts_idx_t *idx)
{
    int actual_n_lvls = idx->n_lvls + 1; /* semantic correction */
    unsigned n_bins;
    int *bins = win2bins(ti, ti + 1, actual_n_lvls, &n_bins);
    unsigned bgzf_bytes = 0;

    /* add initial chunks for [ti, ti + 1) */
    int i;
    for (i = 0; i != n_bins; ++i)
        bgzf_bytes += 
            add_bin_chunks_aux(bins[i], idx->bidx[tid], min_off, itree);
    free(bins);
    return bgzf_bytes;
}


/* transfer all chunks in itree to chunks, appending from initial
   position *n_chunks and updating both chunks and n_chunks. */
static void
accumulate_chunks(kbtree_t(itree_t) *itree,
                  hts_pair64_t **chunks,
                  unsigned *n_chunks)
{
    unsigned beg = *n_chunks;
    *n_chunks += kb_size(itree);
    *chunks = realloc(*chunks, *n_chunks * sizeof((*chunks)[0]));

    kbitr_t tree_itr;
    kb_itr_first_itree_t(itree, &tree_itr);
    int iret = 1; /* */
    unsigned i;
    for (i = beg; i != *n_chunks; ++i) {
        assert(iret != 0);
        (*chunks)[i] = kb_itr_key(hts_pair64_t, &tree_itr);
        iret = kb_itr_next_itree_t(itree, &tree_itr);
    }
}

struct virt_less_range_par {
    struct contig_region *ary;
    struct contig_pos q;
};

/* uses par both as a source of elements elem and a query element q.
   return 1 if q < elem[pos], 0 otherwise.
 */
static int
less_rng_end(unsigned pos, void *par)
{
    struct virt_less_range_par *vl = par;
    struct contig_pos end_pos = { vl->ary[pos].tid, vl->ary[pos].end };
    return cmp_contig_pos(vl->q, end_pos) < 0;
}

/* using forward scanning through the index, return the largest
   position end such that the total size of bgzf blocks overlapping
   [beg, end) is at or slightly above min_wanted_bytes. use [qbeg,
   qend) to define the subset of ranges to consider.  max_end defines
   the maximum logical end range to consider.

   (*chunks)[c] = chunk c
   (*n_chunks) = number of chunks
*/
static struct contig_pos
hts_size_to_range(struct contig_pos beg,
                  struct contig_pos max_end,
                  unsigned min_wanted_bytes,
                  struct contig_region *qbeg,
                  struct contig_region *qend,
                  hts_idx_t *idx,
                  hts_pair64_t **chunks,
                  unsigned *n_chunks,
                  unsigned *more_input)
{
    /* constants */
    int firstbin0 = hts_bin_first(idx->n_lvls);
    int firstbin1 = hts_bin_first(idx->n_lvls - 1);
    int n_small_bins = firstbin0 - firstbin1;
    int n_bins, *bins = malloc(idx->n_lvls * sizeof(int));

    /* initialize qcur */
    struct virt_less_range_par vpar = { qbeg, beg };
    struct contig_region *qcur = qbeg
        + virtual_upper_bound(0, qend - qbeg, less_rng_end, &vpar);
        
    int ms = idx->min_shift; /* convert between 16kb window index and
                                position */
    struct contig_pos cpos; /* current contig position */
    struct contig_region creg;
    unsigned bgzf_bytes;

    for (*n_chunks = bgzf_bytes = 0; qcur != qend; ++qcur) {
        creg = (struct contig_region){ 
            qcur->tid, MAX(qcur->beg, beg.pos), qcur->end
        };
        int ti, b;

        /* creg overlaps positions [creg.beg, creg.end), and overlaps
           16kb windows [ti_beg, ti_end).  if creg.end == 16384, we
           want ti_end == 1.  but if creg.end == 16385, we want ti_end
           = 2. (creg.end - 1) is the last position that creg
           contains, and thus (creg.end - 1)>>ms is the 16kb window
           index of the window containing this position.  we thus add
           1 to this last containing window to get the bound. */
        int ti_beg = creg.beg>>ms,
            ti_end = ((creg.end - 1)>>ms) + 1;
        uint64_t min_off = find_min_offset(creg.tid, ti_beg, idx);
        kbtree_t(itree_t) *itree = kb_init(itree_t, 128);
        bgzf_bytes += compute_initial_bins(creg.tid, ti_beg, min_off, itree, idx);

        cpos.tid = creg.tid;

        /* we use ti = ti_beg + 1 because compute_initial_bins gives
           the chunks for tiling windows [ti, ti + 1). */
        for (ti = ti_beg + 1; ti != ti_end; ++ti) {
            cpos.pos = (ti + 1)<<ms; /* boundary position (one greater
                                        than last position covered by
                                        ti)*/
            if (bgzf_bytes >= min_wanted_bytes
                || cmp_contig_pos(cpos, max_end) == 1) {
                /* break if window touches max_end, or we have enough
                   content */
                accumulate_chunks(itree, chunks, n_chunks);
                kb_destroy(itree_t, itree);
                goto END;
            }
            
            /* process this 16kb window */
            n_bins = next_bins_aux(ti, bins, firstbin0, n_small_bins);
            for (b = 0; b != n_bins; ++b)
                bgzf_bytes += 
                    add_bin_chunks_aux(bins[b], idx->bidx[creg.tid], min_off, itree);
        }
        accumulate_chunks(itree, chunks, n_chunks);
        kb_destroy(itree_t, itree);

        /* beg is only relevant within the first value of qcur */
        beg = (struct contig_pos){ 0, 0 }; 
    }
 END:
    *more_input = qcur != qend;
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
    struct bam_scanner_info *bp = par;
    float mul = STEP_DOWN;
    unsigned tmp_bytes, min_actual_bytes = UINT_MAX, min_wanted_bytes;
    struct contig_pos tmp_end, min_end = { UINT_MAX, UINT_MAX };
    do {
        min_wanted_bytes = max_bytes * mul;
        unsigned s;
        for (s = 0; s != bp->n; ++s) {
            struct bam_stats *bs = &bp->m[s];
            tmp_end =
                hts_size_to_range(cs_stats.pos,
                                  min_end,
                                  min_wanted_bytes,
                                  bp->qbeg,
                                  bp->qend,
                                  bs->idx, 
                                  &bs->chunks,
                                  &bs->n_chunks,
                                  &bs->more_input);
            min_end = MIN_CONTIG_POS(min_end, tmp_end);
            tmp_bytes = tally_bgzf_bytes(bs->chunks, bs->n_chunks);
            min_actual_bytes = MIN(min_actual_bytes, tmp_bytes);
        }
        mul *= STEP_DOWN;
    } while (min_actual_bytes > max_bytes);
    bp->loaded_span = (struct contig_span){ cs_stats.pos, min_end };
    cs_stats.pos = min_end;
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
   calls return the correct number of characters read.  return 1 if
   more input is still available, 0 if this is the last of it. */
unsigned
bam_reader(void *par, struct managed_buf *bufs)
{
    struct bam_scanner_info *bp = par;
    unsigned s, c;
    int block_span, bgzf_block_len, remain;
    ssize_t n_chars_read;
    char *wp;

    unsigned max_grow;
    unsigned more_input = 0;

    for (s = 0; s != bp->n; ++s) {
        struct bam_stats *bs = &bp->m[s];
        if (bs->more_input) more_input = 1;

        wp = bufs[s].buf + bufs[s].size;
        max_grow = tally_bgzf_bytes(bs->chunks, bs->n_chunks);
        ALLOC_GROW_REMAP(bufs[s].buf, wp, bufs[s].size + max_grow, bufs[s].alloc);

        hFILE *fp = bs->bgzf->fp;
        for (c = 0; c != bs->n_chunks; ++c) {
            bgzf_seek(bs->bgzf, bs->chunks[c].u, SEEK_SET);

            /* read all but possibly the last block */
            block_span = BAM_COFFSET(bs->chunks[c].v) - BAM_COFFSET(bs->chunks[c].u);
            n_chars_read = hread(fp, wp, block_span);
            wp += block_span;
            
            if (BAM_UOFFSET(bs->chunks[c].v) != 0) {
                /* the chunk end virtual offset refers to a non-zero
                   offset into a block. */

                /* determine size of last block and read remainder */
                n_chars_read = hread(fp, wp, BLOCK_HEADER_LENGTH);
                
                bgzf_block_len = bgzf_block_size(wp);
                wp += BLOCK_HEADER_LENGTH;
                remain = bgzf_block_len - BLOCK_HEADER_LENGTH;
                
                n_chars_read = hread(fp, wp, remain);
                assert(n_chars_read == remain);
                
                wp += remain;
            }
        }
        bufs[s].size = wp - bufs[s].buf;
    }
    return more_input;
}

static size_t
inflate_bgzf_block(char *in, int block_length, char *out);


/* bgzf contains the bgzf blocks to be inflated.  stores inflated BAM
   records in bam.  manage the size of bam.  only copy the portions of
   blocks defined by chunks and n_chunks. */
void
bam_inflate(const struct managed_buf *bgzf,
            hts_pair64_t *chunks,
            unsigned n_chunks,
            struct managed_buf *bam)
{
    unsigned c;
    int64_t c_beg, c_end; /* sub-regions in a chunk that we want */
    size_t n_copy;
    uint64_t bs; /* compressed block size */

    bam->size = 0;
    char *in = bgzf->buf, *out = bam->buf;
    size_t n_bytes;

/* current and end coffset as we traverse the chunks within our block
   ranges */
    uint64_t coff, coff_beg, coff_end, uoff_end; 
    for (c = 0; c != n_chunks; ++c) {
        coff_beg = BAM_COFFSET(chunks[c].u);
        coff_end = BAM_COFFSET(chunks[c].v); /* offset of last block, if it exists  */
        uoff_end = BAM_UOFFSET(chunks[c].v); /* if non-zero, block at coff_end is present */
        coff = coff_beg;
        /* inflate all blocks starting in [coff_beg, coff_end).
           optionally, inflate block starting at coff_end if it has
           non-zero uoffset (i.e. there are records we need from it) */
        while (coff < coff_end || (coff == coff_end && uoff_end != 0)) {
            ALLOC_GROW_REMAP(bam->buf, out, bam->size + BGZF_MAX_BLOCK_SIZE, bam->alloc);
            bs = bgzf_block_size(in);
            
            n_bytes = inflate_bgzf_block(in, bs, out);
            c_beg = (coff == coff_beg ? BAM_UOFFSET(chunks[c].u) : 0);
            c_end = (coff == coff_end ? uoff_end : n_bytes);
            n_copy = c_end - c_beg; /* number of bytes to copy to output buffer */

            if (c_beg != 0) memmove(out, out + c_beg, n_copy);

            /* test the first record in this block */
            bam1_t brec;
            if (coff == coff_beg) bam_parse(out, &brec);
            // assert(brec.core.tid == 0);

            bam->size += n_copy;
            in += bs;
            out += n_copy;
            coff += bs;
        }
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
    uint32_t *x = (uint32_t *)(raw + 4);
#endif

    c->tid = x[0]; c->pos = x[1];
    c->bin = x[2]>>16; c->qual = x[2]>>8&0xff; c->l_qname = x[2]&0xff;
    c->flag = x[3]>>16; c->n_cigar = x[3]&0xffff;
    c->l_qseq = x[4];
    c->mtid = x[5]; c->mpos = x[6]; c->isize = x[7];

    b->l_data = b->m_data = block_len - 32;
    return raw + block_len + 4;
}


/* allocate a strndup'ed buffer of raw's contents of just one BAM
   record. */
char *
bam_duplicate_buf(char *raw)
{
    /* size of BAM record */
    int32_t n = (*(int32_t *)raw) + 4;
    /* char *dup = strndup(raw, n); */
    char *dup = malloc(n);
    memcpy(dup, raw, n);
    return dup;
    /* !!! for some bizarre reason, using strndup instead of malloc /
           memcpy doesn't seem to work. */
}

/* initialize bs from bam_file.  also, expects <bam_file>.bai to exist
   as well */
void
bam_stats_init(const char *bam_file, struct bam_stats *bs)
{
    bs->bgzf = bgzf_open(bam_file, "r");
    if (! bs->bgzf) {
        fprintf(stderr, "Error: Couldn't load bam file %s\n", bam_file);
        exit(1);
    }
    char *bam_index_file = malloc(strlen(bam_file) + 5);
    strcpy(bam_index_file, bam_file);
    strcat(bam_index_file, ".bai");

    bs->idx = bam_index_load(bam_index_file);
    if (! bs->idx) {
        fprintf(stderr, "Error: Couldn't load bam index file %s\n", bam_index_file);
        exit(1);
    }
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
