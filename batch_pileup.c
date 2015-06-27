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


/* hash type for storing bqt (basecall, quality, strand) counts at
   position */
KHASH_MAP_INIT_INT64(pbqt_h, unsigned);

/* hash type for storing basecall counts at position */
KHASH_MAP_INIT_INT(p_h, unsigned[4]);

/* hash type for storing distinct indels */
KHASH_MAP_INIT_INT(indel_h, struct indel_seq *);

/* hash type for tallying indel events */
KHASH_MAP_INIT_INT(indel_ct_h, khint_t);


__thread struct tls {
    khash_t(indel_h) *indel_hash;
    khash_t(pbqt_h) **pbqt_hash;
    khash_t(p_h) **p_hash;
    khash_t(indel_ct_h) **indel_ct_hash;
    unsigned *tally_beg; /* beginning of range of positions still
                            being tallied. tally_beg[s] for sample
                            s. */
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
    khash_t(pbqt_h) *h_src = tls.pbqt_hash[s];
    khash_t(p_h) *h_trg = tls.p_hash[s];
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


struct pos_base_count {
    unsigned pos;
    unsigned basecount[4];
};

int pos_base_count_less(struct pos_base_count a, struct pos_base_count b)
{
    return a.pos < b.pos;
}

KSORT_INIT(pbc_sort, struct pos_base_count, pos_base_count_less);

/* create a sorted array from sample s's p_hash */
struct pos_base_count *make_p_array(unsigned s)
{
    khash_t(p_h) *ph = tls.p_hash[s];
    unsigned n = kh_size(ph);
    struct pos_base_count *ary = malloc(n * sizeof(struct pos_base_count));
    khint_t k;
    unsigned i;
    for (k = kh_begin(ph), i = 0; k != kh_end(ph); ++k)
        if (kh_exist(ph, k))
            ary[i++] = (struct pos_base_count){ kh_key(ph, k), kh_val(ph, k) };
    
    ks_introsort(pbc_sort, n, ary);
    return ary;
}

struct pos_indel_count {
    unsigned pos;
    khint_t indel_key;
    unsigned count;
};

/* create a sorted array from sample s's indel_hash */
struct pos_indel_count *make_indel_array(unsigned s)
{
    khash_t 
    struct pos_indel_count *ary;
}

