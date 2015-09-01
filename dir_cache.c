/* Utilities for caching bounds and points */

#include "dir_cache.h"
#include "khash.h"
#include "thread_queue.h"
#include "bam_reader.h"
#include "bam_sample_info.h"
#include "dir_diff_cache.h"
#include "batch_pileup.h"
#include "chunk_strategy.h"

#include <pthread.h>
#include <assert.h>

/* Will hold all cached points */
static POINT *g_point_sets_buf;

/* use union alpha_large_key as key */
KHASH_MAP_INIT_INT64(counts_h, unsigned);
static khash_t(counts_h) *g_ptup_hash;
static khash_t(counts_h) *g_btup_hash;

/* protect global hashes */
static pthread_mutex_t merge_mtx = PTHREAD_MUTEX_INITIALIZER;

/* use {pack,unpack}_alpha64 to convert khint64_t key */
KHASH_MAP_INIT_INT64(points_h, POINT *);
khash_t(points_h) *g_points_hash;

/* use {pack,unpack}_bounds to convert khint64_t key */
KHASH_MAP_INIT_INT64(bounds_h, struct binomial_est_bounds);
khash_t(bounds_h) *g_bounds_hash;

/* use {pack,unpack}_alpha64 to convert khint64_t key.  the key
   represents the alpha counts of the sample, and its fuzzy_state
   relationship with the alpha counts { pseudo_depth, 0, 0, 0 } */
KHASH_MAP_INIT_INT64(ref_change_h, enum fuzzy_state);
khash_t(ref_change_h) *g_ref_change_hash;

                                                                                          

static struct dir_cache_params g_dc_par;


#define MIN(a, b) ((a) < (b) ? (a) : (b))
#define MAX(a, b) ((a) < (b) ? (b) : (a))

void
dir_cache_init(struct dir_cache_params dc_par)
{
    g_dc_par = dc_par;
    g_ptup_hash = kh_init(counts_h);
    g_btup_hash = kh_init(counts_h);
    g_points_hash = kh_init(points_h);
    g_bounds_hash = kh_init(bounds_h);
    g_ref_change_hash = kh_init(ref_change_h);
    size_t buf_sz = 
        (size_t)g_dc_par.n_point_sets
        * (size_t)g_dc_par.max_sample_points
        * sizeof(POINT);

    g_point_sets_buf = malloc(buf_sz);
    if (! g_point_sets_buf) {
        fprintf(stderr, "%s:%u: Couldn't allocate point sets buffer of %Zu bytes\n",
                __FILE__, __LINE__, buf_sz);
        exit(1);
    }
}


void
dir_cache_free()
{
    kh_destroy(counts_h, g_ptup_hash);
    kh_destroy(counts_h, g_btup_hash);
    kh_destroy(points_h, g_points_hash);
    kh_destroy(bounds_h, g_bounds_hash);
    kh_destroy(ref_change_h, g_ref_change_hash);
    free(g_point_sets_buf);
}



/* component sizes (24, 20, 12, 8) in bits */
union pair {
    uint32_t c[2];
    uint64_t v;
};


uint64_t
pack_alpha64(unsigned a0, unsigned a1, unsigned a2, unsigned a3)
{
    union pair p;
    p.c[0] = (uint32_t)a0<<8 | (uint32_t)a3;
    p.c[1] = (uint32_t)a1<<12 | (uint32_t)a2;
    return p.v;
}


void
unpack_alpha64(uint64_t k, unsigned *c)
{
    union pair p;
    p.v = k;
    c[0] = p.c[0]>>8;
    c[3] = p.c[0] & 0xff;
    c[1] = p.c[1]>>12;
    c[2] = p.c[1] & 0xfff;
}



static const unsigned alpha_packed_limits[] = {
    1<<24, 1<<20, 1<<12, 1<<8
};

/* defines the limits for using pack_bounds. */
static const unsigned bounds_packed_limits[] = {
    1<<20, 1<<24, 1<<20
};


/* pack 20|24|20 bits of a2|b1|b2.  */
uint64_t
pack_bounds(unsigned a2, unsigned b1, unsigned b2)
{
    uint64_t k = (uint64_t)a2<<44 | (uint64_t)b1<<20 | (uint64_t)b2;
    return k;
}


unsigned
try_pack_bounds(unsigned a2, unsigned b1, unsigned b2, uint64_t *key)
{
    unsigned packable = 
        (a2 <= bounds_packed_limits[0]
         && b1 <= bounds_packed_limits[1]
         && b2 <= bounds_packed_limits[2]);

    if (packable)
        *key = pack_bounds(a2, b1, b2);
    return packable;
}


/* */
void
unpack_bounds(uint64_t k, unsigned *b)
{
    b[0] = k>>44;
    b[1] = (unsigned)(k>>20 & (uint64_t)0xffffff);
    b[2] = (unsigned)(k & (uint64_t)0x0fffff);
}


/* attempts to initialize the key from the counts. returns 1 on
   success, 0 on failure. if counts do not fit within the packing
   limits, does not modify key. */
static unsigned 
try_pack_alpha_aux(unsigned *cts, uint64_t *key)
{
    unsigned packable = 
        cts[0] <= alpha_packed_limits[0]
        && cts[1] <= alpha_packed_limits[1]
        && cts[2] <= alpha_packed_limits[2]
        && cts[3] <= alpha_packed_limits[3];

    if (packable)
        *key = pack_alpha64(cts[0], cts[1], cts[2], cts[3]);
    return packable;
}


/* Given a[4], a_lim[4], b[4], and b_lim[4], find a permutation of
   {(a[p_1], b[p_1]), (a[p_2], b[p_2]), (a[p_3], b[p_3]), (a[p_4],
   b[p_4]) } such that both a[p_i] and b[p_i] are less than their
   corresponding lim[p_i], if such permutation exists.
   
   Assumes that lim[i] >= lim[i+1].  Return 1 if permutation found, 0
   otherwise.  */
static unsigned
find_cacheable_perm_aux(const unsigned *a, const unsigned *a_lim, 
                        const unsigned *b, const unsigned *b_lim,
                        unsigned *permutation)
{
    unsigned perm_found = 1;
    /* mpi[i] (max permutation index) is the maximum position in the
       permutation (described above) that the (a, b) pair may attain
       and still stay below the limits. */
    /* min value, min_index */
    int i, j;
    int min_mpi, mpi[] = { -1, -1, -1, -1 };
    
    unsigned tot[] = { a[0] + b[0], a[1] + b[1], a[2] + b[2], a[3] + b[3] };
    
    for (i = 0; i != 4; ++i)
        for (j = 0; j != 4; ++j) {
            if (a[i] >= a_lim[j] || b[i] >= b_lim[j]) break;
            mpi[i] = j;
        }
    
    /* p[i], (the permutation of mpi), s.t. mpi[p[i]] <= mpi[p[i+1]] */

    /* mpi[j] = -1 means the value j has been placed in the
       permutation */
    for (i = 0; i != 4; ++i) {
        unsigned max_v = 0, max_i = 4;
        min_mpi = 5;
        for (j = 0; j != 4; ++j)
            if (mpi[j] >= i && mpi[j] <= min_mpi && tot[j] >= max_v) {
                max_v = tot[j];
                max_i = j;
                min_mpi = mpi[j];
            }
        if (max_i != 4)
            permutation[i] = max_i, mpi[max_i] = -1;
        else {
            perm_found = 0;
            break;
        }
    }
    return perm_found;
}


/* find the permutation of the sample-to-REF alpha pair. if a
   permutation is found, it is safe to apply pack_alpha64 or
   pack_bounds on the permuted 'a' and/or 'b'. */
unsigned
sample_to_ref_perm(const unsigned *a, const unsigned *b, unsigned *perm)
{
    unsigned pseudo_limits[] = { g_dd_par.pseudo_depth, 1, 1, 1 };
    return 
        find_cacheable_perm_aux(a, alpha_packed_limits,
                                b, pseudo_limits,
                                perm);
}


/* find the permutation of the sample-to-sample alpha pair. if a
   permutation is found, it is safe to apply pack_alpha64 or
   pack_bounds on the permuted 'a' and/or 'b'. */
unsigned
sample_to_sample_perm(const unsigned *a, const unsigned *b, unsigned *perm)
{
    return 
        find_cacheable_perm_aux(a, alpha_packed_limits,
                                b, alpha_packed_limits,
                                perm);
}


POINT *
dir_cache_try_get_points(unsigned *alpha)
{
    khint64_t key = pack_alpha64(alpha[0], alpha[1], alpha[2], alpha[3]);
    khiter_t itr = kh_get(points_h, g_points_hash, key);
    return 
        (itr != kh_end(g_points_hash) && kh_exist(g_points_hash, itr))
        ? kh_val(g_points_hash, itr)
        : NULL;
}


struct binomial_est_bounds *
dir_cache_try_get_bounds(unsigned a2, unsigned b1, unsigned b2)
{
    khint64_t key = pack_bounds(a2, b1, b2);
    khiter_t itr;
    if ((itr = kh_get(bounds_h, g_bounds_hash, key)) != kh_end(g_bounds_hash)
        && kh_exist(g_bounds_hash, itr))
        return &kh_val(g_bounds_hash, itr);
    else return NULL;
}




enum fuzzy_state *
dir_cache_try_get_refchange(unsigned *alpha)
{
    khint64_t key = pack_alpha64(alpha[0], alpha[1], alpha[2], alpha[3]);
    khiter_t itr = kh_get(ref_change_h, g_ref_change_hash, key);
    if (itr != kh_end(g_ref_change_hash) && kh_exist(g_ref_change_hash, itr))
        return &kh_val(g_ref_change_hash, itr);
    else return NULL;
}


/* set the fuzzy state from the cache and return 1.  or, return 0 if
   no cache entry is found. */
unsigned
dir_cache_try_get_diff(unsigned *a, unsigned *b, enum fuzzy_state *state)
{
    unsigned i, c, perm[4], pa[4];
    c = sample_to_ref_perm(a, b, perm);
    if (c) {
        for (i = 0; i != 4; ++i) pa[i] = a[perm[i]];
        *state = *dir_cache_try_get_refchange(pa);
        return 1;
    }
    c = sample_to_sample_perm(a, b, perm);
    if (c) {
        struct binomial_est_bounds *beb =
            dir_cache_try_get_bounds(a[perm[1]], b[perm[0]], b[perm[1]]);
        if (beb) {
            *state = locate_cell(beb, a[perm[0]]);
            return 1;
        }
    }
    return 0;
}




static void
survey_worker(const struct managed_buf *in_bufs,
              unsigned more_input,
              void *scan_info,
              struct managed_buf *out_bufs);

static void
winnow_ptup_hash();


static void
survey_on_create()
{
    batch_pileup_thread_init(bam_samples.n,
                             g_dc_par.fasta_file);
}


static void
survey_on_exit()
{
    batch_pileup_thread_free();
}


/* a no-op.  (all work is done at exit) */
static void
survey_offload(void *par,
               const struct managed_buf *bufs)
{
}




/* main routine for running the survey.  read chunks of input as
   determined by reader_pars until either n_max_loci are read, or at
   least n_bounds and n_point_sets. */
void
run_survey(struct bam_filter_params bf_par,
           struct bam_scanner_info *reader_buf,
           unsigned pseudo_depth,
           unsigned long n_max_survey_loci,
           unsigned n_threads,
           unsigned n_max_reading,
           unsigned long max_input_mem)
{
    unsigned skip_empty_loci = 1;
    batch_pileup_init(bf_par, skip_empty_loci, pseudo_depth);

    thread_queue_reader_t reader = { bam_reader, bam_scanner };
    unsigned t, nb = 0, np = 0;
    khiter_t itr;

    struct contig_span target_span;
    target_span.beg = CONTIG_REGION_BEG(cs_stats.query_regions[0]);
    unsigned long n_loci_left = n_max_survey_loci / bam_samples.n;
    
#define N_LOCI_PER_CHUNK 1e7

    /* need at least min_n_loci in order to keep all threads
       occupied. 16384 is the minimum loci parseable from a BAM
       file. */
    unsigned min_n_loci = 16384 * n_threads;
    unsigned max_n_loci_loop = MAX(min_n_loci, N_LOCI_PER_CHUNK / bam_samples.n);
    unsigned long n_found_loci;
    void **reader_pars = malloc(n_threads * sizeof(void *));
    for (t = 0; t != n_threads; ++t)
        reader_pars[t] = &reader_buf[t];

    /* loop until we have enough statistics or run out of input.  nb,
     np, n_loci_left, target_span, and qcur are all updated in the
     iteration. */
    while ((nb < g_dc_par.n_bounds || np < g_dc_par.n_point_sets)
           && n_loci_left > 0) {
        max_n_loci_loop = MIN(max_n_loci_loop, n_loci_left);
        target_span.end = 
            find_span_of_size(cs_stats.query_regions,
                              cs_stats.query_regions + cs_stats.n_query_regions,
                              target_span.beg, max_n_loci_loop, &n_found_loci);

        if (cmp_contig_pos(target_span.beg, target_span.end) == 0) break;
        
        n_loci_left -= n_found_loci;
        chunk_strategy_set_span(target_span);
            
        struct thread_queue *tq =
            thread_queue_init(reader,
                              reader_pars,
                              survey_worker,
                              survey_offload, NULL,
                              survey_on_create,
                              survey_on_exit,
                              n_threads,
                              0, /* don't need extra output buffers */
                              n_max_reading,
                              bam_samples.n,
                              0, /* no outputs */
                              max_input_mem);
        
        thread_queue_run(tq);
        thread_queue_free(tq);

        target_span.beg = target_span.end;

        nb = 0;
        for (itr = kh_begin(g_btup_hash); itr != kh_end(g_btup_hash); ++itr) {
            if (! kh_exist(g_btup_hash, itr)) continue;
            if (kh_val(g_btup_hash, itr) >= g_dc_par.min_ct_keep_bound) ++nb;
        }
        np = kh_size(g_ptup_hash);
    }
    free(reader_pars);
    batch_pileup_free();
    winnow_ptup_hash();
}



/* merge */
static void
survey_merge(khash_t(counts_h) *ptup_hash,
             khash_t(counts_h) *btup_hash)
{
    khiter_t itr1, itr2;
    khint64_t key;
    unsigned ct;
    int ret;
    for (itr1 = kh_begin(ptup_hash); itr1 != kh_end(ptup_hash); ++itr1) {
        if (! kh_exist(ptup_hash, itr1)) continue;
        key = kh_key(ptup_hash, itr1);
        ct = kh_val(ptup_hash, itr1);
        itr2 = kh_put(counts_h, g_ptup_hash, key, &ret);
        if (ret == 0)
            kh_val(g_ptup_hash, itr2) += ct;
        else
            kh_val(g_ptup_hash, itr2) = ct;
    }
    for (itr1 = kh_begin(btup_hash); itr1 != kh_end(btup_hash); ++itr1) {
        if (! kh_exist(btup_hash, itr1)) continue;
        key = kh_key(btup_hash, itr1);

        ct = kh_val(btup_hash, itr1);
        itr2 = kh_put(counts_h, g_btup_hash, key, &ret);
        if (ret == 0)
            kh_val(g_btup_hash, itr2) += ct;
        else
            kh_val(g_btup_hash, itr2) = ct;
    }
}


#define MIN(a, b) ((a) < (b) ? (a) : (b))

/* delete the lowest-count tuples in g_ptup_hash until the total size is
   less than g_dc_par.n_point_sets. (not a very efficient method -- it
   does multiple passes through the hash unnecessarily) */
static void
winnow_ptup_hash()
{
    unsigned ct, cur_min = 1, nxt_min = UINT_MAX;
    khiter_t itr;
    while (kh_size(g_ptup_hash)) {
        for (itr = kh_begin(g_ptup_hash); itr != kh_end(g_ptup_hash); ++itr) {
            if (kh_size(g_ptup_hash) <= g_dc_par.n_point_sets) return;
            if (! kh_exist(g_ptup_hash, itr)) continue;
            ct = kh_val(g_ptup_hash, itr);
            if (ct == cur_min) kh_del(counts_h, g_ptup_hash, itr);
            else nxt_min = MIN(ct, nxt_min);
        }
        cur_min = nxt_min;
        nxt_min = UINT_MAX;
    }
}


/* accumulate basecall statistics */
static void
survey_worker(const struct managed_buf *in_bufs,
              unsigned more_input,
              void *scan_info,
              struct managed_buf *out_bufs /* unused */)
{
    struct managed_buf bam = { NULL, 0, 0 };
    struct bam_scanner_info *bsi = scan_info;
    unsigned s;
    pileup_load_refseq_ranges(bsi);
    for (s = 0; s != bam_samples.n; ++s) {
        struct bam_stats *bs = &bsi->m[s];
        bam_inflate(&in_bufs[s], bs->chunks, bs->n_chunks, &bam);
        pileup_tally_stats(bam, bsi, s);
    }
    if (bam.buf != NULL) free(bam.buf);

    if (! more_input)
        pileup_final_input();

    /* we need bqs and indel stats.  do not need basecall stats */
    for (s = 0; s != bam_samples.n; ++s)
        pileup_prepare_basecalls(s);

    struct locus_data 
        *ld[2], 
        *ldat = malloc(bam_samples.n * sizeof(struct locus_data)),
        pseudo_data;

    for (s = 0; s != bam_samples.n; ++s)
        alloc_locus_data(&ldat[s]);

    alloc_locus_data(&pseudo_data);

    khint64_t key;
    unsigned perm[4], perm_found;
    khiter_t itr;
    int ret;
    unsigned i, pi;
    unsigned *ct[2];
    khash_t(counts_h) *ptup_hash = kh_init(counts_h);
    khash_t(counts_h) *btup_hash = kh_init(counts_h);

    while (pileup_next_pos()) {
        reset_locus_data(&pseudo_data);
        for (s = 0; s != bam_samples.n; ++s)
            reset_locus_data(&ldat[s]);

        for (pi = 0; pi != bam_sample_pairs.n; ++pi) {
            unsigned sp[] = { bam_sample_pairs.m[pi].s1, 
                              bam_sample_pairs.m[pi].s2 };

            ld[0] = &ldat[sp[0]];
            ld[1] = sp[1] == PSEUDO_SAMPLE ? &pseudo_data : &ldat[sp[1]];

            for (i = 0; i != 2; ++i)
                if (! ld[i]->init.base_ct) {
                    ld[i]->base_ct = pileup_current_basecalls(sp[i]);
                    ld[i]->init.base_ct = 1;
                }

            ct[0] = ld[0]->base_ct.ct_filt;
            ct[1] = ld[1]->base_ct.ct_filt;

            if (sp[1] != PSEUDO_SAMPLE) {
                perm_found = sample_to_sample_perm(ct[0], ct[1], perm);
                if (perm_found) {
                    /* tally the bounds tuple */
                    key = pack_bounds(ct[0][perm[1]], ct[1][perm[0]], ct[1][perm[1]]);
                    itr = kh_put(counts_h, btup_hash, key, &ret);
                    if (ret == 0) kh_val(btup_hash, itr)++;
                    else kh_val(btup_hash, itr) = 1;
                }
            }

            /* store individual sample alphas. */
            if (perm_found)
                for (i = 0; i != 2; ++i) {
                    key = pack_alpha64(ct[i][perm[0]], ct[i][perm[1]], 
                                       ct[i][perm[2]], ct[i][perm[3]]);
                    itr = kh_put(counts_h, ptup_hash, key, &ret);
                    if (ret == 0) kh_val(ptup_hash, itr)++;
                    else kh_val(ptup_hash, itr) = 1;
                }

        }

    }
    pileup_clear_stats();

    for (s = 0; s != bam_samples.n; ++s)
        free_locus_data(&ldat[s]);
    free(ldat);
    free_locus_data(&pseudo_data);
    
    pthread_mutex_lock(&merge_mtx);
    survey_merge(ptup_hash, btup_hash);
    pthread_mutex_unlock(&merge_mtx);

    kh_destroy(counts_h, ptup_hash);
    kh_destroy(counts_h, btup_hash);

    /* fprintf(stdout, "Surveyed %u:%u-%u\t%u\n",  */
    /*         bsi->loaded_span.beg.tid, bsi->loaded_span.beg.pos, */
    /*         bsi->loaded_span.end.pos, */
    /*         bsi->loaded_span.end.pos - bsi->loaded_span.beg.pos); */
}


struct gen_points_par {
    unsigned n_threads, nth;
};


static void *
generate_points_worker(void *par)
{
    struct gen_points_par gp = *(struct gen_points_par *)par;
    struct points_gen_par pgp = {
        .alpha_perm = { 0, 1, 2, 3 },
        .randgen = gsl_rng_alloc(gsl_rng_taus)
    };
    
    khiter_t itr, itr2;
    unsigned i;
    khint64_t key;
    POINT *points, *pend;
    size_t block_ct = g_dc_par.max_sample_points;
    int ret;

    khash_t(points_h) *ph = kh_init(points_h);
    for (i = 0, itr = kh_begin(g_ptup_hash); itr != kh_end(g_ptup_hash); ++itr) {
        if (! kh_exist(g_ptup_hash, itr)) continue;
        if (i % gp.n_threads == gp.nth) {
            /* create the hash entry */
            key = kh_key(g_ptup_hash, itr);
            points = g_point_sets_buf + (i * block_ct);
            itr2 = kh_put(points_h, ph, key, &ret);
            assert(ret != 0);
            kh_val(ph, itr2) = points;

            /* prepare for points generation */
            unpack_alpha64(key, pgp.alpha_counts);

            /* generate the points */
            pend = points + block_ct;
            while (points != pend) {
                gen_dirichlet_points_wrapper(&pgp, points);
                points += GEN_POINTS_BATCH;
            }
        }
        ++i;
    }
    gsl_rng_free(pgp.randgen);

    /* merge into final points hash */
    pthread_mutex_lock(&merge_mtx);
    for (itr = kh_begin(ph); itr != kh_end(ph); ++itr) {
        if (! kh_exist(ph, itr)) continue;
        key = kh_key(ph, itr);
        points = kh_val(ph, itr);
        itr2 = kh_put(points_h, g_points_hash, key, &ret);
        assert(ret != 0);
        kh_val(g_points_hash, itr2) = points;
    }
    pthread_mutex_unlock(&merge_mtx);
    kh_destroy(points_h, ph);
    return NULL;
}



/* Populate g_points_hash with the total set of point sets from the
   dirichlets of g_ptup_hash */
void
generate_point_sets(unsigned n_threads)
{
    pthread_t *threads = malloc(n_threads * sizeof(pthread_t));
    struct gen_points_par *gpp = malloc(n_threads * sizeof(struct gen_points_par));
    unsigned t;
    for (t = 0; t != n_threads; ++t) {
        gpp[t] = (struct gen_points_par){ n_threads, t };
        pthread_create(&threads[t], NULL, generate_points_worker, &gpp[t]);
    }
    for (t = 0; t != n_threads; ++t)
        pthread_join(threads[t], NULL);

    free(gpp);
    free(threads);
}


struct gen_est_bounds_par {
    unsigned n_threads, nth;
};

static void *
generate_est_bounds_worker(void *par)
{
    struct gen_est_bounds_par geb = *(struct gen_est_bounds_par *)par;
    khiter_t itr, itr2;
    khash_t(bounds_h) *bh = kh_init(bounds_h);
    khint64_t key;

    unsigned i, bounds[3];
    int ret;
    struct binomial_est_bounds beb;
    
    struct bound_search_params bpar = {
        .dist = { malloc(sizeof(struct distrib_points)),
                  malloc(sizeof(struct distrib_points)) }
    };

    alloc_distrib_points(bpar.dist[0]);
    alloc_distrib_points(bpar.dist[1]);

    for (i = 0, itr = kh_begin(g_btup_hash); itr != kh_end(g_btup_hash); ++itr) {
        if (! kh_exist(g_btup_hash, itr)) continue;
        if (i % geb.n_threads == geb.nth) {
            /* create key */
            key = kh_key(g_btup_hash, itr);
            itr2 = kh_put(bounds_h, bh, key, &ret);
            assert(ret != 0);
            unpack_bounds(key, bounds);
            initialize_est_bounds(bounds[0], bounds[1], bounds[2],
                                  &bpar, &beb);

            kh_val(bh, itr2) = beb;
        }
        ++i;
    }

    free_distrib_points(bpar.dist[0]);
    free_distrib_points(bpar.dist[1]);
    free(bpar.dist[0]);
    free(bpar.dist[1]);

    /* add to global hash */
    pthread_mutex_lock(&merge_mtx);
    for (itr = kh_begin(bh); itr != kh_end(bh); ++itr) {
        if (! kh_exist(bh, itr)) continue;
        key = kh_key(bh, itr);
        beb = kh_val(bh, itr);
        itr2 = kh_put(bounds_h, g_bounds_hash, key, &ret);
        assert(ret != 0);
        kh_val(g_bounds_hash, itr2) = beb;
    }
    
    pthread_mutex_unlock(&merge_mtx);
    kh_destroy(bounds_h, bh);
    return NULL;
}


/* populate g_bounds_hash with computed est bounds for all of the
   bounds tuples surveyed. */
void
generate_est_bounds(unsigned n_threads)
{
    pthread_t *threads = malloc(n_threads * sizeof(pthread_t));
    struct gen_est_bounds_par *geb = 
        malloc(n_threads * sizeof(struct gen_est_bounds_par));

    unsigned t;
    for (t = 0; t != n_threads; ++t) {
        geb[t] = (struct gen_est_bounds_par){ n_threads, t };
        pthread_create(&threads[t], NULL, generate_est_bounds_worker, &geb[t]);
    }
    for (t = 0; t != n_threads; ++t)
        pthread_join(threads[t], NULL);

    free(geb);
    free(threads);
}


/* populate g_ref_change_hash with computed fuzzy_states */
