/* Provides ability to populate a rolling buffer worth of pileup
   statistics for multiple samples in tandem, operating on one BAM
   record at a time and updating all single locus statistics that it
   affects.

Approach:

"Tally phase"
for each sample: 
   for each BAM record:
       update endpos[s] to start of BAM record
       for each position:
           add to PBQT hash
           add to main indel hash
           add to PK hash

   Create a PB hash from the PBQT hash
   Create a PB array from PB hash
   Sort PB array by position

   Create a PK array from the PK hash
   Sort PK array by position

"Summary Phase"
for each sample pair:
    for each position in pair of pb_ct arrays:
        compute the cached difference
            if above-threshold
                compute weighted difference using PBQT hash (enumerate all
                possible keys at position)


this is the clear_finished() function:

Delete all entries in each PBQT hash that fall before the next
processed position.

Delete all entries in PK hashes that fall before next processed
position

Delete PB arrays
Delete PK arrays
Clear PB hashes



NOTE: The residual entries left over in each PBQT hash and PK hash
contain the information from the BAM records for positions that are
not yet complete.  By saving this information, we can avoid doing
re-tally on any BAM records.  They will become complete once more BAM
records are processed.

During the tally phase, each sample's PBQT hash takes on average 160
bytes per position, and is populated fully all in one tight loop,
allowing good cache usage.  The PK hash will be much smaller but
follows the same approach.

Then, the two hashes are efficiently traversed from start to finish
and summarized into more compact structures.  The resulting P hash has
~1000 positions.  The array compiled from it then is sorted, which,
for such a small array, is negligible cost.

More importantly, the array is now only about 5kb, and pairs of these
arrays will very efficiently fit in the smallest caches during the
summary phase, where they are used together in pairs.

thread-local storage will be used for these hashes for convenience.

*/

struct indel_seq {
    char is_ins;
    char seq[FLEX_ARRAY];
};

struct pbqt {
    uint32_t pos;
    unsigned : 22;
    unsigned long base: 2;
    unsigned long qual: 7;
    unsigned long strand: 1;
};

union pbqt_key {
    khint64_t k;
    struct pbqt v;
};

struct pos_base_count {
    unsigned pos;
    unsigned basecount[4];
};


struct indel_count {
    khint_t key;
    unsigned ct;
};

struct pos_indel_count {
    unsigned pos;
    struct indel_count i;
};

/* hash type for storing bqt (basecall, quality, strand) counts at
   position */
KHASH_MAP_INIT_INT64(pbqt_h, unsigned);

/* hash type for storing basecall counts at position */
KHASH_MAP_INIT_INT(p_h, unsigned[4]);

/* hash type for storing distinct indels at position */
KHASH_MAP_INIT_INT(indel_h, struct indel_seq *);

/* hash type for tallying indel events */
KHASH_MAP_INIT_INT(indel_ct_h, struct indel_count);


struct tally_stats {
    khash_t(pbqt_h) *pbqt_hash;
    khash_t(indel_ct_h) *indel_ct_hash;

    khash_t(p_h) *p_hash;
    struct pos_base_count *base_ct_ary;
    struct pos_indel_count *indel_ct_ary;

    unsigned tally_beg; /* beginning of range of positions still being
                           tallied. */
};

__thread struct tls {
    khash_t(indel_h) *indel_hash;
    unsigned n_samples;
    struct tally_stats *ts; /* ts[s] is the tally stats for a given sample */
    unsigned cur_summary_pos; /* current position being summarized */
} tls;


/* tally one BAM record for sample s, updating tls hashes. return
   pointer to next BAM record. */
char *tally_bam(unsigned s, char *bam)
{
    /* traverse the BAM record */
    
}


/* condense to p_hash */
void summarize_base_counts(unsigned s)
{
    union pbqt_key k;
    unsigned ct;
    khash_t(pbqt_h) *h_src = tls.ts[s].pbqt_hash;
    khash_t(p_h) *h_trg = tls.ts[s].p_hash;
    kh_clear(p_h, h_trg);
    khint_t i, j;
    for (i = kh_begin(h_src); i != kh_end(h_src); ++i)
    {
        if (! kh_exist(h_src, i)) continue;
        k.k = kh_key(h_src, i);
        ct = kh_val(h_src, i);
        j = kh_get(h_trg, k.v.pos);
        if (j == kh_end(h_trg)) 
        {
            j = kh_put(p_h, h_trg, k.v.pos);
            kh_val(h_trg, j) = { 0, 0, 0, 0 };
        }
        kh_val(h_trg, j)[k.base] += ct;
    }
}


int pos_base_count_less(struct pos_base_count a, struct pos_base_count b)
{
    return a.pos < b.pos;
}

KSORT_INIT(pbc_sort, struct pos_base_count, pos_base_count_less);

/* create a sorted array from sample s's p_hash */
void make_p_array(unsigned s)
{
    khash_t(p_h) *ph = tls.ts[s].p_hash;
    unsigned n = kh_size(ph);
    struct pos_base_count *ary = malloc(n * sizeof(struct pos_base_count));

    khint_t k;
    unsigned i;
    for (k = kh_begin(ph), i = 0; k != kh_end(ph); ++k)
        if (kh_exist(ph, k))
            ary[i++] = (struct pos_base_count){ kh_key(ph, k), kh_val(ph, k) };
    
    ks_introsort(pbc_sort, n, ary);
    tls.ts[s].base_ct_ary = ary;
}

int pos_indel_count_less(struct pos_indel_count a, struct pos_indel_count b)
{
    return a.pos < b.pos;
}

KSORT_INIT(pi_sort, struct pos_indel_count, pos_indel_count_less);

/* create a sorted array from sample s's indel_hash */
void make_indel_array(unsigned s)
{
    khash_t(indel_h) *ih = tls.ts[s].indel_hash;
    unsigned n = kh_size(ih);
    struct pos_indel_count *ary = malloc(n * sizeof(struct pos_indel_count));
    khint_t k;
    unsigned i;
    for (k = kh_begin(ih), i = 0; k != kh_end(ih); ++k)
        if (kh_exist(ih, k))
            ary[i++] = (struct pos_indel_count){ kh_key(ih, k), kh_val(ih, k) };

    ks_introsort(pi_sort, n, ary);
    tls.ts[s].indel_ct_ary = ary;
}


/* clear statistics that are no longer needed */
void clear_finished_stats()
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
            if (pk.v.pos < tls.cur_summary_pos) kh_del(h, i);
        }

        /* clear finished entries in indel counts hash */
        khash_t(indel_ct_h) *ih = tls.ts[s].indel_ct_hash;
        for (i = kh_begin(ih); i != kh_end(ih); ++i)
        {
            if (! kh_exist(ih, i)) continue;
            pos = kh_key(ih, i);
            if (pos < tls.cur_summary_pos) kh_del(ih, i);
        }

        /* completely clear temporary data */
        kh_clear(p_h, tls.ts[s].p_hash);
        free(tls.ts[s].base_ct_ary);
        free(tls.ts[s].indel_ct_ary);
    }
}


/* sample pairwise calculations */
