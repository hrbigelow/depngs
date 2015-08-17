/* Utilities for caching bounds and points */

#include "khash.h"
#include "thread_queue.h"


/* use union alpha_large_key as key */
KHASH_MAP_INIT_INT64(points_tuple_h, unsigned);

/* use union bounds_key as key */
KHASH_MAP_INIT_INT64(bounds_tuple_h, unsigned);

/* use union alpha_large_key as key */
KHASH_MAP_INIT_INT64(points_h, double *);

/* use union bounds_key as key */
KHASH_MAP_INIT_INT64(bounds_h, struct binomial_est_bounds);

/* protect global hashes */
static pthread_mutex_t merge_mtx = PTHREAD_MUTEX_INITIALIZER;

/* Survey */
static khash_t(points_tuple_h) *g_pt_hash;
static khash_t(bounds_tuple_h) *g_bt_hash;

/* Prepopulation */
static khash_t(points_h) *g_points_hash;
static khash_t(bounds_h) *g_bounds_hash;

/* Will hold all cached points */
static double *g_point_sets_buf;
                                                                                          
{};



/* main routine for running the survey.  read chunks of input until
   accumulating at least n_bounds and n_point_sets.  At each
   iteration,  */
void
run_survey(void **reader_pars,
           struct contig_region *qbeg,
           struct contig_region *qend,
           unsigned n_threads,
           unsigned n_max_reading,
           unsigned long max_input_mem,
           unsigned n_bounds,
           unsigned min_ct_keep_bound,
           unsigned n_point_sets)
{
    thread_queue_reader_t reader = { bam_reader, bam_scanner };
    unsigned nb = 0, np = 0;
    khiter_t itr;
    unsigned t;

    /* clone the original logical range */
    unsigned n_ranges = qend - qbeg;
    assert(n_ranges);

    struct contig_region
        *regs = malloc(n_ranges * sizeof(struct contig_region)),
        *breg = regs,
        *ereg = regs + n_ranges,
        *creg;
    memcpy(regs, qbeg, n_ranges);
    unsigned save_end = breg->end;

    /* parse one million data points per chunk. A data point is one
       locus in one sample */
#define N_DATA_PER_CHUNK 1e6

    unsigned max_n_loci = N_DATA_PER_CHUNK / bam_samples.n;
    struct bam_scanner_info *bsi;

    while (nb < n_bounds || np < n_point_sets) {

        /* scan forward until we have n_loci */
        unsigned n_loci = 0;
        breg->beg = breg->end;
        breg->end = save_end;
        for (breg = creg; breg != ereg; ++breg) {
            n_loci += breg->end - breg->beg;
            if (n_loci > max_n_loci) {
                save_end = breg->end;
                breg->end -= (n_loci - max_n_loci);
            }
        }
            
        /* set the next chunk of work for each thread */
        for (t = 0; t != n_threads; ++t) {
            bsi = reader_pars[t];
            bsi->qbeg = breg;
            bsi->qend = ereg;
        }
    
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

        nb = 0;
        for (itr = kh_begin(g_bounds_hash); itr != kh_end(g_bounds_hash); ++itr) {
            if (! kh_exist(g_bounds_hash, itr)) continue;
            if (kh_val(g_bounds_hash, itr) >= min_ct_keep_bound) ++nb;
        }
        np = kh_size(g_points_hash);
    }
    free(regs);

    /* restore original settings for the reader_pars */
    for (t = 0; t != n_threads; ++t) {
        bsi = reader_pars[t];
        bsi->qbeg = qbeg;
        bsi->qend = qend;
    }

}


static void
survey_on_create()
{
}


static void
survey_on_exit()
{
}


/* a no-op.  (all work is done at exit) */
static void
survey_offload(void *par,
               const struct managed_buf *bufs)
{
}


/* merge */
static void
survey_merge(khash_t(points_tuple_h) *pt_hash,
             khash_t(bounds_tuple_h) *bt_hash)
{
    khiter_t itr1, itr2;
    union alpha_large_key ak;
    unsigned ct;
    for (itr1 = kh_begin(pt_hash); itr1 != kh_end(pt_hash); ++itr1) {
        if (! kh_exist(pt_hash, itr1)) continue;
        ak.raw = kh_key(pt_hash, itr1);
        ct = kh_val(pt_hash, itr1);
        itr2 = kh_put(points_tuple_h, g_pt_hash, ak.raw, &ret);
        if (ret == 0)
            kh_val(g_pt_hash, itr2) += ct;
        else
            kh_val(g_pt_hash, itr2) = ct;
    }
    union bounds_key bk;
    for (itr1 = kh_begin(bt_hash); itr1 != kh_end(bt_hash); ++itr1) {
        if (! kh_exist(bt_hash, itr1)) continue;
        bk.val = kh_key(bt_hash, itr1);
        ct = kh_val(bt_hash, itr1);
        itr2 = kh_put(bounds_tuple_h, g_bt_hash, bk.val, &ret);
        if (ret == 0)
            kh_val(g_bt_hash, itr2) += ct;
        else
            kh_val(g_bt_hash, itr2) = ct;
    }
}


/* accumulate basecall statistics */
static void
survey_worker(const struct managed_buf *in_bufs,
              unsigned more_input,
              void *scan_info,
              struct managed_buf *out_bufs)
{
    struct managed_buf bam = { NULL, 0, 0 };
    struct bam_scanner_info *bsi = vsi;
    unsigned s;
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

    struct locus_data ld[2];
    alloc_locus_data(&ld[0]);
    alloc_locus_data(&ld[1]);

    union alpha_large_key ak;
    union bounds_key bk;
    unsigned perm[4], perm_found, *cts;
    int ret;
    unsigned i, pi;
    khash_t(points_tuple_h) *pt_hash = kh_init(points_tuple_h);
    khash_t(bounds_tuple_h) *bt_hash = kh_init(bounds_tuple_h);

    while (pileup_next_pos()) {
        for (pi = 0; pi != bam_sample_pairs.n; ++pi) {
            unsigned sp[] = { bam_sample_pairs.m[pi].s1, 
                              bam_sample_pairs.m[pi].s2 };

            for (i = 0; i != 2; ++i) {
                if (! ld[i]->init.base_ct) {
                    ld[i]->base_ct = pileup_current_basecalls(sp[i]);
                    ld[i]->init.base_ct = 1;
                }
            }
            /* tally the dirichlet base count */
            find_cacheable_permutation(ld[0]->base_ct.ct_filt,
                                       ld[1]->base_ct.ct_filt,
                                       alpha_packed_limits,
                                       perm,
                                       &perm_found);
            if (! perm_found) continue;
            
            for (i = 0; i != 2; ++i) {
                cts = ld[i]->base_ct.ct_filt;
                ak.c = (struct alpha_packed_large){ 
                    cts[perm[0]],
                    cts[perm[1]],
                    cts[perm[2]],
                    cts[perm[3]]
                }
                itr = kh_put(points_tuple_h, pt_hash, ak.raw, &ret);
                if (ret == 0)
                    kh_val(pt_hash, itr)++;
                else
                    kh_val(pt_hash, itr) = 1;
            }

            /* tally the bounds tuple */
            bk.f.a2 = ld[0]->base_ct.ct_filt[perm[1]];
            bk.f.b1 = ld[1]->base_ct.ct_filt[perm[0]];
            bk.f.b2 = ld[1]->base_ct.ct_filt[perm[1]];

            itr = kh_put(bounds_tuple_h, bt_hash, bk.val, &ret);
            if (ret == 0)
                kh_val(bt_hash, itr)++;
            else
                kh_val(bt_hash, itr) = 1;
        }

        free_locus_data(&ld[0]);
        free_locus_data(&ld[1]);
    }

    pthread_mutex_lock(&merge_mtx);
    survey_merge(pt_hash, bt_hash);
    pthread_mutex_unlock(&merge_mtx);

    kh_destroy(pt_hash);
    kh_destroy(bt_hash);
}


struct gen_points_par {
    unsigned n_threads, nth;
    unsigned max_sample_points;
};

static void *
generate_points_worker(void *par)
{
    struct gen_points_par gp = *(struct gen_points_par *)par;
    struct points_gen_par pgp = {
        .alpha_perm = { 0, 1, 2, 3 },
        .randgen = gsl_rng_alloc(gsl_rng_taus);
    };

    
    struct points_gen pgen = { &pgp, gen_dirichlet_points_wrapper, NULL };
    
    khiter_t itr, itr2;
    unsigned i;
    union alpha_large_key ak;
    double *points, *pend;
    size_t block_ct = 4 * gp.max_sample_points;

    khash_t(points_h) *ph = kh_init(points_h);
    for (i = 0, itr = kh_begin(g_pt_hash); itr != kh_end(g_pt_hash); ++itr) {
        if (! kh_exist(g_pt_hash, itr)) continue;
        if (i % gp.n_threads == gp.nth) {
            /* create the hash entry */
            ak.raw = kh_key(g_pt_hash, itr);
            points = g_point_sets_buf + (i * block_ct);
            itr2 = kh_put(points_h, ph, ak.raw, &ret);
            assert(t != 0);
            kh_val(ph, itr2) = points;

            /* prepare for points generation */
            pgp.alpha_counts[0] = ak.c.a0;
            pgp.alpha_counts[1] = ak.c.a1;
            pgp.alpha_counts[2] = ak.c.a2;
            pgp.alpha_counts[3] = ak.c.a3;

            /* generate the points */
            pend = points + block_ct;
            while (points != pend) {
                gen_dirichlet_points_wrapper(&pgp, points);
                points += GEN_POINTS_BATCH * 4;
            }
        }
        ++i;
    }
    gsl_rng_free(pgp.randgen);

    /* merge into final points hash */
    pthread_mutex_lock(&merge_mtx);
    for (itr = kh_begin(ph); itr != kh_end(ph); ++itr) {
        if (! kh_exist(ph, itr)) continue;
        ak.raw = kh_key(ph, itr);
        points = kh_val(ph, itr);
        itr2 = kh_put(points_h, g_points_hash, ak.raw, &ret);
        assert(ret != 0);
        kh_val(g_points_hash, itr2) = points;
    }
    pthread_mutex_unlock(&merge_mtx);
    kh_destroy(ph);
}



/* Populate g_points_hash with the total set of point sets from the
   dirichlets of g_pt_hash */
void
generate_point_sets(unsigned n_threads,
                    unsigned max_sample_points)
{
    unsigned t;
    pthread_t *threads = malloc(n_threads * sizeof(pthread_t));
    struct gen_points_par gpp = { n_threads, 0, max_sample_points };
    for (t = 0; t != n_threads; ++t) {
        gpp.nth = t;
        pthread_create(&threads[t], NULL, generate_points_worker, &gpp);
    }
    for (t = 0; t != n_threads; ++t)
        pthread_join(&threads[t], NULL);

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
    union bounds_key bk;

    unsigned i;
    struct binomial_est_bounds beb;
    struct binomial_est_params bpar;

    for (i = 0, itr = kh_begin(g_bt_hash); itr != kh_end(g_bt_hash); ++itr) {
        if (! kh_exist(g_bt_hash, itr)) continue;
        if (i % geb.n_threads == geb.nth) {
            /* create key */
            bk.val = kh_key(g_bt_hash, itr);
            itr2 = kh_put(bounds_h, bh, bk.key, &ret);
            assert(ret != 0);
            
            initialize_est_bounds(bk.f.a2, bk.f.b1, bk.f.b2,
                                  &bpar, &beb);

            kh_val(bh, itr2) = beb;
        }
        ++i;
    }

    /* add to global hash */
    pthread_mutex_lock(&merge_mtx);
    for (itr = kh_begin(bh); itr != kh_end(bh); ++itr) {
        if (! kh_exist(bh, itr)) continue;
        bk.val = kh_key(bh, itr);
        beb = kh_val(bh, itr);
        itr2 = kh_put(bounds_h, g_bounds_hash, bk.val, &ret);
        assert(ret != 0);
        kh_val(g_bounds_hash, itr2) = beb;
    }
    
    pthread_mutex_unlock(&merge_mtx);
    kh_destroy(bh);
}


/* populate g_bounds_hash with computed est bounds for all of the
   bounds tuples surveyed. */
void
generate_est_bounds(unsigned n_threads)
{
    unsigned t;
    pthread_t *threads = malloc(n_threads * sizeof(pthread_t));
    struct gen_est_bounds_par geb = { n_threads, 0 };

    for (t = 0; t != n_threads; ++t) {
        geb.nth = t;
        pthread_create(&threads[t], NULL, generate_est_bounds_worker, &geb);
    }
    for (t = 0; t != n_threads; ++t)
        pthread_join(&threads[t], NULL);

    free(threads);
}
