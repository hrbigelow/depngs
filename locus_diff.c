#define __STDC_LIMIT_MACROS
#include <stdint.h>

#include "locus_diff.h"

#include "sampling.h"

#include "yepLibrary.h"
#include "yepCore.h"

#include <gsl/gsl_math.h>
#include <gsl/gsl_randist.h>
#include <time.h>
#include <pthread.h>
#include <assert.h>

#include "bam_reader.h"
#include "batch_pileup.h"
#include "khash.h"
#include "bam_sample_info.h"
#include "common_tools.h"
#include "locus.h"
#include "chunk_strategy.h"
#include "geometry.h"
#include "dirichlet_points_gen.h"


#define MIN(a, b) ((a) < (b) ? (a) : (b))

static struct timespec start_time;
static double posterior_confidence;
static double min_dirichlet_dist;
static unsigned max_sample_points;
static double indel_alpha_prior;

static struct {
    unsigned do_print_pileup;
    unsigned do_dist, do_comp, do_indel;
} worker_options;


static struct {
    struct bam_scanner_info *reader_buf;
    void **reader_pars;
    unsigned n_max_reading;
    struct contig_region *ranges;
    unsigned n_ranges;
    unsigned n_threads;
    struct locus_diff_offload_par offload_par;
    const char *fasta_file;
} thread_params;


static __thread struct locus_diff_input tls_dw;


static double quantiles[MAX_NUM_QUANTILES];
static unsigned n_quantiles;


void
alloc_locus_data(struct locus_data *ld)
{
    alloc_distrib_points(&ld->distp);
    ld->bqs_ct = NULL;
    ld->indel_ct = NULL;
    init_pileup_data(&ld->sample_data);
}


void
free_locus_data(struct locus_data *ld)
{
    free_distrib_points(&ld->distp);
    free(ld->bqs_ct);
    free(ld->indel_ct);
    free_pileup_data(&ld->sample_data);
}


/* call when we advance to a new locus */
void
reset_locus_data(struct locus_data *ld)
{
    ld->init.distp = 0;
    ld->init.base_ct = 0;
    ld->init.bqs_ct = 0;
    ld->init.indel_ct = 0;
    ld->init.sample_data = 0;
    ld->distp.points.size = 0;
    ld->distp.weights.size = 0;
    ld->confirmed_changed = 0;
}


void
locus_diff_init(double _post_confidence, 
                double _beta_confidence,
                double _min_dirichlet_dist,
                unsigned _max_sample_points,
                unsigned _max_dir_cache_items,
                unsigned _max_bounds_cache_items,
                unsigned n_threads,
                double prior_alpha,
                const char *samples_file,
                const char *sample_pairs_file,
                const char *fasta_file,
                unsigned min_quality_score,
                const char *quantiles_string,
                unsigned do_dist,
                unsigned do_comp,
                unsigned do_indel,
                unsigned do_print_pileup)
{
    clock_gettime(CLOCK_REALTIME, &start_time);
    if (_post_confidence < 0.8 || _post_confidence > 0.999999) {
        fprintf(stderr,
                "Error: locus_diff_init: posterior confidence of %g is a bad value\n",
                _post_confidence);
        exit(1);
    }
    posterior_confidence = _post_confidence;
    if (_min_dirichlet_dist <= 0 || _min_dirichlet_dist >= 1) {
        fprintf(stderr,
                "Error: locus_diff_init: min_dirichlet_dist of %g is a bad value.\n"
                "Should be in [0, 1]\n", _min_dirichlet_dist);
        exit(1);
    }
    min_dirichlet_dist = _min_dirichlet_dist;
    max_sample_points = _max_sample_points;

    bam_sample_info_init(samples_file, sample_pairs_file);

    /* we do not want to skip empty loci, because we need to traverse
       these in order to get statistics for missing data */
    unsigned skip_empty_loci = 0;
    batch_pileup_init(min_quality_score, skip_empty_loci);

    dirichlet_diff_cache_init(PSEUDO_DEPTH,
                              GEN_POINTS_BATCH,
                              _post_confidence, 
                              _beta_confidence,
                              prior_alpha,
                              _min_dirichlet_dist,
                              _max_sample_points,
                              _max_dir_cache_items, 
                              _max_bounds_cache_items,
                              n_threads);

    parse_csv_line(quantiles_string, quantiles, &n_quantiles, MAX_NUM_QUANTILES);

    worker_options.do_print_pileup = do_print_pileup;
    worker_options.do_dist = do_dist;
    worker_options.do_comp = do_comp;
    worker_options.do_indel = do_indel;
}


void
locus_diff_free()
{
    bam_sample_info_free();
    batch_pileup_free();
    dirichlet_diff_cache_free();
}


/* called just after this thread is created */
void
dist_on_create()
{
    tls_dw.randgen = gsl_rng_alloc(gsl_rng_taus);
    alloc_locus_data(&tls_dw.pseudo_sample);

    tls_dw.ldat = malloc(bam_samples.n * sizeof(struct locus_data));
    unsigned s;
    for (s = 0; s != bam_samples.n; ++s)
        alloc_locus_data(&tls_dw.ldat[s]);

    tls_dw.pair_stats = calloc(bam_sample_pairs.n, sizeof(struct pair_dist_stats));
    tls_dw.square_dist_buf = malloc(sizeof(double) * max_sample_points);
    tls_dw.weights_buf = malloc(sizeof(double) * max_sample_points);
    tls_dw.do_print_progress = 1; /* !!! how to choose which thread prints progress? */
    tls_dw.bep.points_hash_frozen = 0;

    batch_pileup_thread_init(bam_samples.n, thread_params.fasta_file);
}


/* when a worker thread exits, it must inform the rest of the program
   that it will not be modifying shared data anymore. */
void
dist_on_exit()
{
    gsl_rng_free(tls_dw.randgen);

    free_locus_data(&tls_dw.pseudo_sample);
    unsigned s;
    for (s = 0; s != bam_samples.n; ++s)
        free_locus_data(&tls_dw.ldat[s]);
    free(tls_dw.ldat);

    free(tls_dw.pair_stats);
    free(tls_dw.square_dist_buf);
    free(tls_dw.weights_buf);

    batch_pileup_thread_free();

    inactivate_shared_data(! tls_dw.bep.points_hash_frozen,
                           ! tls_dw.bep.bounds_hash_frozen);

}


#define NUM_EXTRA_BUFS_PER_THREAD 500

struct thread_queue *
locus_diff_tq_init(const char *locus_range_file,
                   const char *fasta_file,
                   unsigned n_threads,
                   unsigned n_max_reading,
                   unsigned long max_input_mem,
                   FILE *dist_fh,
                   FILE *comp_fh,
                   FILE *indel_fh)
{
    thread_params.n_threads = n_threads;

    unsigned long n_total_loci;
    thread_params.ranges = 
        parse_locus_ranges(locus_range_file,
                           fasta_file,
                           &thread_params.n_ranges,
                           &n_total_loci);

    thread_params.reader_buf = 
        malloc(n_threads * sizeof(struct bam_scanner_info));

    thread_params.reader_pars = malloc(n_threads * sizeof(void *));

    unsigned t, s;
    for (t = 0; t != n_threads; ++t) {
        thread_params.reader_buf[t] = (struct bam_scanner_info){
            malloc(bam_samples.n * sizeof(struct bam_stats)),
            bam_samples.n,
            
            thread_params.ranges, 
            thread_params.ranges + thread_params.n_ranges
        };
        for (s = 0; s != bam_samples.n; ++s)
            bam_stats_init(bam_samples.m[s].bam_file, 
                           &thread_params.reader_buf[t].m[s]);
        
        thread_params.reader_pars[t] = &thread_params.reader_buf[t];
    }
    thread_params.fasta_file = fasta_file;

    if (locus_range_file)
        cs_init_by_range(n_total_loci, bam_samples.n);

    else {
        cs_init_whole_file(bam_samples.n);
        /* for (s = 0; s != bam_samples.n; ++s) */
        /*     cs_set_total_bytes(s, ix[s].root->end_offset - ix[s].root->start_offset); */
    }

#define MAX_BYTES_SMALL_CHUNK 1e9
#define SMALL_CHUNK 5e6
#define DEFAULT_BYTES_PER_LOCUS 100

    cs_set_defaults(MAX_BYTES_SMALL_CHUNK,
                    SMALL_CHUNK, 
                    DEFAULT_BYTES_PER_LOCUS);

    thread_params.offload_par = 
        (struct locus_diff_offload_par){ dist_fh, comp_fh, indel_fh };

    /* To avoid a stall, n_extra / n_threads should be greater than
       Max(work chunk time) / Avg(work chunk time). */
    size_t n_extra = n_threads * NUM_EXTRA_BUFS_PER_THREAD;
    size_t n_output_files = 
        (dist_fh ? 1 : 0) + (comp_fh ? 1 : 0) + (indel_fh ? 1 : 0);

    thread_queue_reader_t reader = { bam_reader, bam_scanner };

    struct thread_queue *tqueue =
        thread_queue_init(reader, 
                          thread_params.reader_pars,
                          locus_diff_worker,
                          locus_diff_offload, 
                          &thread_params.offload_par,
                          dist_on_create,
                          dist_on_exit,
                          n_threads, n_extra, n_max_reading, bam_samples.n,
                          n_output_files, max_input_mem);

    return tqueue;
}


void
locus_diff_tq_free()
{
    unsigned t, s;
    for (t = 0; t != thread_params.n_threads; ++t) {
        for (s = 0; bam_samples.n; ++s)
            bam_stats_free(&thread_params.reader_buf[t].m[s]);

        free(thread_params.reader_buf[t].m);
    }

    free(thread_params.reader_buf);
    free(thread_params.reader_pars);
    free(thread_params.ranges);

}

typedef double COMP_QV[NUM_NUCS][MAX_NUM_QUANTILES];

struct dim_mean {
    unsigned dim;
    double mean;
};

int
less_mean(const void *ap, const void *bp)
{
    const struct dim_mean *a = ap, *b = bp;
    return a->mean < b->mean;
}
 
void
print_basecomp_quantiles(COMP_QV quantile_values,
                         const double *means,
                         size_t n_quantiles,
                         const char *label_string,
                         struct pileup_locus_info *pli,
                         struct pileup_data *pdat,
                         struct managed_buf *mb)
{
    char line_label[2048];

    sprintf(line_label,
            "%s\t%s\t%i\t%c\t%u\t%u\t%u",
            label_string, 
            pli->refname,
            pli->pos,
            pli->refbase,
            pdat->n_match_hi_q,
            pdat->n_match_lo_q,
            pdat->n_indel);

    struct dim_mean dim_to_mean[NUM_NUCS];

    unsigned d;
    for (d = 0; d != NUM_NUCS; ++d)
        dim_to_mean[d] = (struct dim_mean){ means[d], d };

    qsort(dim_to_mean, NUM_NUCS, sizeof(dim_to_mean[0]), less_mean);

    /* calculate mean rank order */
    size_t mean_rank_order[NUM_NUCS];
    for (d = 0; d != NUM_NUCS; ++d)
        mean_rank_order[dim_to_mean[d].dim] = d;

    unsigned grow = NUM_NUCS * (sizeof(line_label) + 30 + (11 * MAX_NUM_QUANTILES));
    ALLOC_GROW_TYPED(mb->buf, mb->size + grow, mb->alloc);

    static const char dimension_labels[] = "ACGT";
    for (d = 0; d != NUM_NUCS; ++d) {
        mb->size += sprintf(mb->buf + mb->size,
                            "%s\t%c\t%Zu\t%10.8f", line_label, 
                            dimension_labels[d], mean_rank_order[d], means[d]);
        
        unsigned q;
        for (q = 0; q != n_quantiles; ++q)
            mb->size += sprintf(mb->buf + mb->size, "\t%10.8f", quantile_values[d][q]);

        mb->size += sprintf(mb->buf + mb->size, "\n");
    }
}


/* populates square_dist_buf with squares of euclidean distances
   between points1 and points2 in the barycentric space (R4,
   normalized positive components).  populates weights_buf with
   product of weights1 and weights2. */
void
compute_wsq_dist(const double *points1,
                 const double *weights1,
                 const double *points2,
                 const double *weights2,
                 size_t n_points,
                 double *square_dist_buf,
                 double *weights_buf)
{
    compute_square_dist(points1, points2, n_points, 4, square_dist_buf);
    (void)yepCore_Multiply_V64fV64f_V64f(weights1, weights2, weights_buf, n_points);
}


/* Compute the requested set of distance quantile values from two sets
   of weighted points. */
void
compute_dist_wquantiles(double *square_dist_buf,
                        double *weights_buf,
                        size_t n_points,
                        const double *quantiles,
                        size_t n_quantiles,
                        double *dist_quantile_values)
{
    compute_marginal_wquantiles(square_dist_buf, weights_buf, n_points, 1, 0,
                                quantiles, n_quantiles,
                                dist_quantile_values);

    unsigned q;
    for (q = 0; q != n_quantiles; ++q) 
        dist_quantile_values[q] = sqrt(dist_quantile_values[q]);
}


/* print out distance quantiles. use a pseudo-sample for the second
   one if its pair index is equal to REFERENCE_SAMPLE */
void
print_distance_quantiles(const char *contig,
                         size_t position,
                         char ref_base,
                         size_t pair_index,
                         struct locus_data *ldat,
                         double *dist_quantile_values,
                         struct managed_buf *buf)
{
    unsigned s[] = { bam_sample_pairs.m[pair_index].s1,
                     bam_sample_pairs.m[pair_index].s2 };

    unsigned space = (2 * MAX_LABEL_LEN) + 3 + 100 + (10 * MAX_NUM_QUANTILES);

    ALLOC_GROW_TYPED(buf->buf, buf->size + space, buf->alloc);

    buf->size +=
        sprintf(buf->buf + buf->size, 
                "%s\t%s\t%s\t%Zu\t%c", 
                bam_samples.m[s[0]].label,
                s[1] == REFERENCE_SAMPLE ? "REF" : bam_samples.m[s[1]].label,
                contig,
                position,
                ref_base);
    
    unsigned q;
    for (q = 0; q != n_quantiles; ++q)
        buf->size += sprintf(buf->buf + buf->size, "\t%7.4f",
                             dist_quantile_values[q]);

    if (worker_options.do_print_pileup) {
        unsigned i;
        for (i = 0; i != 2; ++i)
            if (! ldat[s[i]].init.sample_data) {
                ldat[s[i]].init.sample_data = 1;
                pileup_current_data(s[i], &ldat[s[i]].sample_data);
            }
        struct pileup_data *pd1 = &ldat[s[0]].sample_data;
        struct pileup_data *pd2 = &ldat[s[1]].sample_data;
        
        unsigned extra_space = 
            pd1->calls.size + pd1->quals.size
            + pd2->calls.size + pd2->quals.size
            + 50;
        
        ALLOC_GROW_TYPED(buf->buf, buf->size + extra_space, buf->alloc);
        
        buf->size += sprintf(buf->buf + buf->size, 
                             "\t%Zu\t%s\t%s\t%Zu\t%s\t%s",
                             pd1->quals.size,
                             pd1->calls.buf,
                             pd1->quals.buf,
                             pd2->quals.size,
                             pd2->calls.buf,
                             pd2->quals.buf);
    }
    buf->size += sprintf(buf->buf + buf->size, "\n");
}


#define ONE_OVER_SQRT2 0.70710678118654752440

/* for each sample pair, calculate whether the loci differ above the
   given level and confidence. if a pair differs, set the
   confirmed_changed flag for each sample in the pair. if out_buf is
   not NULL, print out distance quantiles.  also generates sample
   points for each sample as needed, both for the preliminary test and
   more points for the final test */
void
distance_quantiles_aux(struct managed_buf *out_buf)
{
    enum fuzzy_state diff_state = AMBIGUOUS;

    unsigned cacheable, cache_was_set;
    struct locus_data *ld[2];
    unsigned pi, i;

    for (pi = 0; pi != bam_sample_pairs.n; ++pi) {
        unsigned sp[] = { bam_sample_pairs.m[pi].s1, bam_sample_pairs.m[pi].s2 };

        ld[0] = &tls_dw.ldat[sp[0]];
        ld[1] = sp[1] == REFERENCE_SAMPLE ? &tls_dw.pseudo_sample : &tls_dw.ldat[sp[1]];
        
        tls_dw.metrics.total++;
        ++tls_dw.pair_stats[pi].total;

        /* load this particular pair of distp into the bep */
        for (i = 0; i != 2; ++i) {
            tls_dw.bep.dist[i] = &ld[i]->distp;
            if (! ld[i]->init.base_ct) {
                ld[i]->base_ct = pileup_current_basecalls(sp[i]);
                ld[i]->init.base_ct = 1;
            }
        }

        diff_state = UNCHANGED;
        diff_state =
            cached_dirichlet_diff(ld[0]->base_ct.ct_filt,
                                  ld[1]->base_ct.ct_filt,
                                  &tls_dw.bep,
                                  &cacheable,
                                  &cache_was_set);

        tls_dw.pair_stats[pi].dist_count[diff_state]++;
        tls_dw.pair_stats[pi].cacheable += cacheable;
        tls_dw.pair_stats[pi].cache_was_set += cache_was_set;

        tls_dw.metrics.cacheable += cacheable;
        tls_dw.metrics.cache_was_set += cache_was_set;

        if (diff_state == CHANGED) {
            /* Finish sampling and do full distance marginal estimation */
            for (i = 0; i != 2; ++i) {
                if (! ld[i]->init.bqs_ct) {
                    ld[i]->init.bqs_ct = 1;
                    pileup_current_bqs(sp[i], &ld[i]->bqs_ct, &ld[i]->n_bqs_ct);
                }

                struct distrib_points *dst = tls_dw.bep.dist[i];
                struct points_gen_par *pgp = dst->pgen.points_gen_par;
                pgp->observed = ld[i]->bqs_ct;
                pgp->n_observed = ld[i]->n_bqs_ct;
                
                /* Generate all points */
                unsigned perm[] = { 0, 1, 2, 3 };
                update_points_gen_params(dst, ld[i]->base_ct.ct_filt, perm);
                POINT *p,
                    *pb = dst->points.buf,
                    *pe = dst->points.buf + max_sample_points;
                for (p = pb; p != pe; p += GEN_POINTS_BATCH)
                    dst->pgen.gen_point(pgp, p);

                dst->points.size = max_sample_points;
                
                /* Generate all weights */
                double *w,
                    *wb = dst->weights.buf,
                    *we = dst->weights.buf + max_sample_points;
                for (w = wb, p = pb; w != we; 
                     w += GEN_POINTS_BATCH, p += GEN_POINTS_BATCH)
                    dst->pgen.weight(p, pgp, w);

                dst->weights.size = max_sample_points;
            }
            /* Compute weighted square distances (max value of 2).
               e.g.
               minimum simplex distance on [0,1.00] scale is 0.25.
               minimum real    distance on [0,1.41] scale is 0.35355
               minimum sq_real distance on [0,2.00] scale is 0.12499
            */

            compute_wsq_dist((const double *)tls_dw.bep.dist[0]->points.buf, 
                             tls_dw.bep.dist[0]->weights.buf,
                             (const double *)tls_dw.bep.dist[1]->points.buf, 
                             tls_dw.bep.dist[1]->weights.buf,
                             max_sample_points,
                             tls_dw.square_dist_buf,
                             tls_dw.weights_buf);

            double test_quantile = 1.0 - posterior_confidence, test_quantile_value;

            /* Compute the test distance quantile (relative to the
               dist, not squared dist) */
            compute_dist_wquantiles(tls_dw.square_dist_buf,
                                    tls_dw.weights_buf,
                                    max_sample_points,
                                    &test_quantile,
                                    1,
                                    &test_quantile_value);

            if (test_quantile_value > min_dirichlet_dist) {
                ++tls_dw.pair_stats[pi].confirmed_changed;
                ld[0]->confirmed_changed = 1;
                ld[1]->confirmed_changed = 1;

                if (out_buf) {
                    compute_dist_wquantiles(tls_dw.square_dist_buf,
                                            tls_dw.weights_buf,
                                            max_sample_points,
                                            quantiles,
                                            n_quantiles,
                                            tls_dw.dist_quantile_values);            
                    
                    /* quantile values are in squared distance terms, and
                       in the [0, 1] x 4 space of points.  The maximum
                       distance between such points is sqrt(2.0).  We want
                       to re-scale it to be 1. */
                    unsigned q;
                    for (q = 0; q != n_quantiles; ++q)
                        tls_dw.dist_quantile_values[q] *= ONE_OVER_SQRT2;
                    
                    struct pileup_locus_info pli;
                    pileup_current_info(&pli);
                    print_distance_quantiles(pli.refname,
                                             pli.pos,
                                             pli.refbase,
                                             pi,
                                             tls_dw.ldat,
                                             tls_dw.dist_quantile_values, 
                                             out_buf);
                }
            }
        }
        
    }
}


void
print_indel_distance_quantiles(size_t pair_index,
                               double *dist_quantile_values,
                               struct indel_pair_count *events,
                               size_t n_events,
                               struct managed_buf *mb)
{
    int s1 = bam_sample_pairs.m[pair_index].s1,
        s2 = bam_sample_pairs.m[pair_index].s2;

    struct locus_data *ld[] = {
        &tls_dw.ldat[s1],
        s2 == REFERENCE_SAMPLE ? &tls_dw.pseudo_sample : &tls_dw.ldat[s2]
    };


    unsigned space = (2 * MAX_LABEL_LEN) + 3 + 100 + (10 * MAX_NUM_QUANTILES);
    ALLOC_GROW_TYPED(mb->buf, mb->size + space, mb->alloc);

    struct pileup_locus_info pli;
    pileup_current_info(&pli);

    mb->size += sprintf(mb->buf + mb->size, 
                        "%s\t%s\t%s\t%c\t%u", 
                        bam_samples.m[s1].label,
                        s2 == REFERENCE_SAMPLE ? "REF" : bam_samples.m[s2].label,
                        pli.refname, pli.refbase, pli.pos);
    
    for (size_t q = 0; q != n_quantiles; ++q)
        mb->size += sprintf(mb->buf + mb->size, "\t%.4f", dist_quantile_values[q]);

    unsigned indel_space = n_events * (10 + 10);
    ALLOC_GROW_TYPED(mb->buf, mb->size + indel_space, mb->alloc);

    /* print indel event counts for each sample */
    struct indel_pair_count *eb = events, *ee = eb + n_events;
    unsigned s;
    for (s = 0, eb = events; s != 2; ++s) {
        mb->size += sprintf(mb->buf + mb->size, "%u",
                            ld[s]->base_ct.n_match_lo_q +
                            ld[s]->base_ct.n_match_hi_q);
        while (eb != ee) {
            mb->size += sprintf(mb->buf + mb->size, ",%i", eb->count[s]);
            eb++;
        }
    }
    eb = events;

    /* retrieve each indel one by one */
    mb->size += sprintf(mb->buf + mb->size, "@");
    static char ins[] = "-+";
    while (eb != ee) {
        struct indel_seq *isq = pileup_current_indel_seq(&eb->indel);
        unsigned grow = 2 + strlen(isq->seq);
        ALLOC_GROW_TYPED(mb->buf, mb->size + grow, mb->alloc);
        mb->size += sprintf(mb->buf + mb->size, ",%c%s", 
                            ins[(unsigned)isq->is_ins], isq->seq);
        eb++;
        free(isq);
    }

    if (worker_options.do_print_pileup) {
        struct pileup_data
            *lp0 = &ld[0]->sample_data,
            *lp1 = &ld[1]->sample_data;

        unsigned extra_space = 
            lp0->calls.size + lp0->quals.size
            + lp1->calls.size + lp1->quals.size
            + 50;

        ALLOC_GROW_TYPED(mb->buf, mb->size + extra_space, mb->alloc);
        mb->size += sprintf(mb->buf + mb->size,
                            "\t%u\t%u\t%u\t%s\t%s\t%u\t%u\t%u\t%s\t%s",
                            lp0->n_match_hi_q,
                            lp0->n_match_lo_q,
                            lp0->n_indel,
                            lp0->calls.buf,
                            lp0->quals.buf,
                            lp1->n_match_hi_q,
                            lp1->n_match_lo_q,
                            lp1->n_indel,
                            lp1->calls.buf,
                            lp1->quals.buf);
    }
    mb->size += sprintf(mb->buf + mb->size, "\n");
}


// print out all next distance quantiles for indels
void
indel_distance_quantiles_aux(struct managed_buf *buf)
{
    struct indel_pair_count *indel_pairs = NULL;

    for (size_t pi = 0; pi != bam_sample_pairs.n; ++pi) {
        int s[] = { bam_sample_pairs.m[pi].s1, bam_sample_pairs.m[pi].s2 };
        
        struct locus_data *ld[] = {
            &tls_dw.ldat[s[0]], 
            s[1] == REFERENCE_SAMPLE ? &tls_dw.pseudo_sample : &tls_dw.ldat[s[1]]
        };
        
        unsigned i;
        for (i = 0; i != 2; ++i) {
            if (! ld[i]->init.indel_ct) {
                ld[i]->init.indel_ct = 1;
                pileup_current_indels(s[i], &ld[i]->indel_ct, &ld[i]->n_indel_ct);
            }
            if (! ld[i]->init.base_ct) {
                ld[i]->init.base_ct = 1;
                ld[i]->base_ct = pileup_current_basecalls(s[i]);
            }
        }
        if (ld[0]->n_indel_ct == 0 && ld[1]->n_indel_ct == 0
            && min_dirichlet_dist > 0) continue;
        
        /* traverse the sorted indel counts in tandem */
        unsigned n_indel_pairs;
        pileup_make_indel_pairs(ld[0]->indel_ct, ld[0]->n_indel_ct,
                                ld[1]->indel_ct, ld[1]->n_indel_ct,
                                &indel_pairs, &n_indel_pairs);
        
        /* need at least some indels */
        if (n_indel_pairs == 0) continue;
        
        /* if # of CIGAR match reads > 0 for either sample, we need
           one extra 'event' */
        unsigned n_match[] = { 
            ld[0]->base_ct.n_match_lo_q + ld[0]->base_ct.n_match_hi_q,
            ld[1]->base_ct.n_match_lo_q + ld[1]->base_ct.n_match_hi_q
        };

        /* in order to compare the two dirichlet posteriors, we
           consider that the CIGAR 'match' event exists if either one
           of the bam_samples has non-zero counts.  If both are zero,
           though, we don't include the CIGAR match event in the
           dirichlet's at all. */
        unsigned has_match_reads = (n_match[0] != 0 || n_match[1] != 0);
        unsigned n_events = n_indel_pairs + has_match_reads;

        double *alpha[] = { malloc(n_events * sizeof(double)),
                            malloc(n_events * sizeof(double))
        };
        size_t buf_sz = n_events * max_sample_points;
        double *points[] = { malloc(buf_sz * sizeof(double)),
                             malloc(buf_sz * sizeof(double))
        };

        unsigned c;
        for (i = 0; i != 2; ++i) {
            for (c = 0; c != n_indel_pairs; ++c)
                alpha[i][c] = indel_pairs[c].count[i] + indel_alpha_prior;
            if (has_match_reads)
                alpha[i][n_indel_pairs] = n_match[i] + indel_alpha_prior;
        }

        gsl_ran_dirichlet(tls_dw.randgen, n_events, alpha[i], points[i]);
        
        compute_square_dist(points[0], points[1], 
                            max_sample_points, n_events,
                            tls_dw.square_dist_buf);

        double test_quantile = 1.0 - posterior_confidence, test_quantile_value;
            
        // compute distance quantiles
        compute_marginal_quantiles(tls_dw.square_dist_buf,
                                   max_sample_points,
                                   1, /* one dimensional */
                                   0, /* use the first dimension */
                                   &test_quantile,
                                   1, /* evaluate only one quantile */
                                   &test_quantile_value);

        test_quantile_value = sqrt(test_quantile_value);
        if (test_quantile_value >= min_dirichlet_dist) {
            compute_marginal_quantiles(tls_dw.square_dist_buf,
                                       max_sample_points,
                                       1, /* one dimensional */
                                       0, /* use the first dimension */
                                       quantiles,
                                       n_quantiles,
                                       tls_dw.dist_quantile_values);
            unsigned q;
            for (q = 0; q != n_quantiles; ++q)
                tls_dw.dist_quantile_values[q] = 
                    sqrt(tls_dw.dist_quantile_values[q]) * ONE_OVER_SQRT2;

            print_indel_distance_quantiles(pi, 
                                           tls_dw.dist_quantile_values,
                                           indel_pairs, n_indel_pairs, buf);
        }

        free(points[0]);
        free(points[1]);
        free(alpha[0]);
        free(alpha[1]);
    }
    free(indel_pairs);
}


void
comp_quantiles_aux(struct managed_buf *comp_buf)
{
    unsigned s;
    struct pileup_locus_info ploc;
    pileup_current_info(&ploc);

    COMP_QV comp_quantile_values;
    double comp_means[NUM_NUCS];

    for (s = 0; s != bam_samples.n; ++s) {
        struct locus_data *ld = &tls_dw.ldat[s];
        
        if (ld->confirmed_changed) {
            unsigned d;
            for (d = 0; d != NUM_NUCS; ++d) {
                compute_marginal_wquantiles((double *)ld->distp.points.buf,
                                            ld->distp.weights.buf,
                                            max_sample_points,
                                            NUM_NUCS,
                                            d,
                                            quantiles,
                                            n_quantiles,
                                            comp_quantile_values[d]);
                comp_means[d] = 
                    compute_marginal_mean((double *)ld->distp.points.buf,
                                          ld->distp.weights.buf,
                                          max_sample_points,
                                          NUM_NUCS,
                                          d);
            }
            if (! ld->init.sample_data) {
                ld->init.sample_data = 1;
                pileup_current_data(s, &ld->sample_data);

                print_basecomp_quantiles(comp_quantile_values,
                                         comp_means,
                                         n_quantiles,
                                         bam_samples.m[s].label,
                                         &ploc,
                                         &ld->sample_data,
                                         comp_buf);
            }
        }
    }
}


/* receives a certain number of in_bufs and a certain number of
   out_bufs.
   
   there is one struct locus_data for each input.  it's current
   field points to the current line being processed. 'gs' is a single
   index indicating the sample with the lowest 'current' among all of
   them.  it is this position that must be fully processed before any
   bam_samples may advance.
*/
void
locus_diff_worker(const struct managed_buf *in_bufs,
                  unsigned more_input,
                  void *vsi,
                  struct managed_buf *out_bufs)
{
    struct timespec worker_start_time;
    clock_gettime(CLOCK_REALTIME, &worker_start_time);
    tls_dw.metrics.total = 0;
    tls_dw.metrics.cacheable = 0;
    tls_dw.metrics.cache_was_set = 0;

    unsigned i = 0;
    struct managed_buf 
        *dist_buf = worker_options.do_dist ? &out_bufs[i++] : NULL,
        *comp_buf = worker_options.do_comp ? &out_bufs[i++] : NULL,
        *indel_buf = worker_options.do_indel ? &out_bufs[i++] : NULL;

    struct managed_buf bam = { NULL, 0, 0 };
    unsigned s;
    struct bam_scanner_info *bsi = vsi;
    for (s = 0; s != bam_samples.n; ++s) {
        struct bam_stats *bs = &bsi->m[s];
        bam_inflate(&in_bufs[s], bs->chunks, bs->n_chunks, &bam);
        pileup_tally_stats(bam, bsi, s);
    }
    if (bam.buf != NULL) free(bam.buf);

    if (! more_input)
        pileup_final_input();

    /* we need all three types of data, so must prepare all of it. */
    for (s = 0; s != bam_samples.n; ++s) {
        pileup_prepare_bqs(s);
        pileup_prepare_basecalls(s);
        pileup_prepare_indels(s);
    }


    /* zero out the pair_stats */
    unsigned pi;
    for (pi = 0; pi != bam_sample_pairs.n; ++pi)
        memset(&tls_dw.pair_stats[pi], 0, sizeof(tls_dw.pair_stats[0]));

    while (pileup_next_pos()) {
        /* this will strain the global mutex, but only during the hash
           loading phase. */
        if (! tls_dw.bep.points_hash_frozen)
            tls_dw.bep.points_hash_frozen = freeze_points_hash();

        if (! tls_dw.bep.bounds_hash_frozen)
            tls_dw.bep.bounds_hash_frozen = freeze_bounds_hash();

        if (dist_buf || comp_buf)
            distance_quantiles_aux(dist_buf);
        
        if (indel_buf)
            indel_distance_quantiles_aux(indel_buf);

        if (comp_buf)
            comp_quantiles_aux(comp_buf);

        for (s = 0; s != bam_samples.n; ++s)
            reset_locus_data(&tls_dw.ldat[s]);
    }   

    /* frees statistics that have already been used in one of the
       distance calculations. */
    pileup_clear_stats();

    if (tls_dw.do_print_progress) {
        struct timespec now;
        clock_gettime(CLOCK_REALTIME, &now);
        unsigned elapsed = now.tv_sec - start_time.tv_sec;

        time_t cal = time(NULL);
        char *ts = strdup(ctime(&cal));
        ts[strlen(ts)-1] = '\0';
        fprintf(stdout, 
                "%s (%02i:%02i:%02i elapsed): Finished processing %s %i\n", 
                ts,
                elapsed / 3600,
                (elapsed % 3600) / 60,
                elapsed % 60,
                "1",
                // fasta_get_contig(bsi->loaded_span.end.tid),
                bsi->loaded_span.end.pos + 1);
        fflush(stdout);
        free(ts);
    }

    accumulate_pair_stats(tls_dw.pair_stats);

    // struct timespec worker_end_time;
    // clock_gettime(CLOCK_REALTIME, &worker_end_time);
    // unsigned elapsed = worker_end_time.tv_sec - worker_start_time.tv_sec;

    // fprintf(stderr, "%Zu\t%u\t%u\t%u\t%u\n", 
    //         tls_dw.thread_num, elapsed, tls_dw.metrics.total,
    //         tls_dw.metrics.cacheable, tls_dw.metrics.cache_was_set);
    
    // print_primary_cache_size();
    // print_cache_stats();
}
 
    
void
locus_diff_offload(void *par, const struct managed_buf *bufs)
{
    struct locus_diff_offload_par *ol = (struct locus_diff_offload_par *)par;
    unsigned i = 0;
    if (ol->dist_fh) {
        fwrite(bufs[i].buf, 1, bufs[i].size, ol->dist_fh), i++;
        fflush(ol->dist_fh);
    }

    if (ol->comp_fh) {
        fwrite(bufs[i].buf, 1, bufs[i].size, ol->comp_fh), i++;
        fflush(ol->comp_fh);
    }

    if (ol->indel_fh) {
        fwrite(bufs[i].buf, 1, bufs[i].size, ol->indel_fh), i++;
        fflush(ol->indel_fh);
    }
}
