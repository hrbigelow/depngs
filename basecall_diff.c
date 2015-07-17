#define __STDC_LIMIT_MACROS
#include <stdint.h>

#include "sampling.h"
#include "yepLibrary.h"
#include "yepCore.h"

#include <gsl/gsl_math.h>
#include <gsl/gsl_randist.h>
#include <time.h>
#include <pthread.h>
#include <assert.h>

#include "basecall_diff.h"
#include "bam_reader.h"
#include "batch_pileup.h"
#include "khash.h"
#include "dist_aux.h"
#include "common_tools.h"
#include "locus.h"
#include "chunk_strategy.h"
#include "geometry.h"
#include "dirichlet_points_gen.h"
#include "indel_diff.h"
#include "pair_stats.h"

#define MIN(a, b) ((a) < (b) ? (a) : (b))

static struct timespec start_time;
static double posterior_confidence;
static double min_dirichlet_dist;
static unsigned max_sample_points;

#define PSEUDO_SAMPLE -1

/* describes all pairs of samples */
static struct sample_pairs sample_pairs;


/* attributes intrinsic to one sample */
static struct sample_attrs sample_attrs;

static struct {
    unsigned do_print_pileup;
    unsigned do_dist, do_comp, do_indel;
} worker_options;


static struct {
    struct bam_reader_par *reader_buf;
    void **reader_par;
    unsigned n_readers;
    struct pair_ordering_range *ranges;
    unsigned n_ranges;
    unsigned n_threads;
    struct dist_worker_offload_par offload_par;
} thread_params;


static __thread struct dist_worker_input tls_dw;


static double quantiles[MAX_NUM_QUANTILES];
static unsigned n_quantiles;


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
}


void
free_locus_data(struct locus_data *ld)
{
    free_distrib_points(&ld->distp);
    free(ld->bqs_ct);
    free(ld->indel_ct);
    free_pileup_data(&ld->sample_data);
}



void dist_worker_init(double _post_confidence, 
                      double _min_dirichlet_dist,
                      unsigned _max_sample_points,
                      const char *samples_file,
                      const char *sample_pairs_file,
                      const char *fastq_type,
                      const char *quantiles_string,
                      unsigned do_dist,
                      unsigned do_comp,
                      unsigned do_indel,
                      unsigned do_print_pileup)
{
    clock_gettime(CLOCK_REALTIME, &start_time);
    if (_post_confidence < 0.8 || _post_confidence > 0.999999)
    {
        fprintf(stderr,
                "Error: dist_worker_init: posterior confidence of %g is a bad value\n",
                _post_confidence);
        exit(1);
    }
    posterior_confidence = _post_confidence;
    if (_min_dirichlet_dist <= 0 || _min_dirichlet_dist >= 1)
    {
        fprintf(stderr,
                "Error: dist_worker_init: min_dirichlet_dist of %g is a bad value.\n"
                "Should be in [0, 1]\n", _min_dirichlet_dist);
        exit(1);
    }
    min_dirichlet_dist = _min_dirichlet_dist;
    max_sample_points = _max_sample_points;

    init_dist_aux();
    sample_pairs = init_sample_pairs(sample_pairs_file);
    sample_attrs = init_sample_attributes(samples_file);
    free_dist_aux();


    parse_csv_line(quantiles_string, quantiles, &n_quantiles);

    worker_options.do_print_pileup = do_print_pileup;
    worker_options.do_dist = do_dist;
    worker_options.do_comp = do_comp;
    worker_options.do_indel = do_indel;
}


void dist_worker_free()
{
    free(sample_pairs.p);
    unsigned s;
    for (s = 0; s != sample_attrs.n; ++s)
    {
        free(sample_attrs.atts[s].file);
        fclose(sample_attrs.atts[s].fh);
    }
    free(sample_attrs.atts);
}


/* called just after this thread is created */
void dist_on_create()
{
    tls_dw.randgen = gsl_rng_alloc(gsl_rng_taus);
    tls_dw.ldat = malloc(sample_attrs.n * sizeof(struct locus_data));

    tls_dw.pair_stats = calloc(sample_pairs.n, sizeof(struct pair_dist_stats));
    tls_dw.square_dist_buf = malloc(sizeof(double) * max_sample_points);
    tls_dw.weights_buf = malloc(sizeof(double) * max_sample_points);
    tls_dw.do_print_progress = 1; /* !!! how to choose which thread prints progress? */
    tls_dw.bep.points_hash_frozen = 0;
}


/* when a worker thread exits, it must inform the rest of the program
   that it will not be modifying shared data anymore. */
void dist_on_exit()
{
    gsl_rng_free(tls_dw.randgen);

    free(tls_dw.ldat);
    free(tls_dw.pair_stats);
    free(tls_dw.square_dist_buf);
    free(tls_dw.weights_buf);

    inactivate_shared_data(! tls_dw.bep.points_hash_frozen,
                           ! tls_dw.bep.bounds_hash_frozen);

}


#define NUM_EXTRA_BUFS_PER_THREAD 500

struct thread_queue *dist_worker_tq_init(const char *query_range_file,
                                         unsigned n_threads,
                                         unsigned n_readers,
                                         unsigned long max_input_mem,
                                         FILE *dist_fh,
                                         FILE *comp_fh,
                                         FILE *indel_fh)
{
    thread_params.n_threads = n_threads;

    thread_params.reader_buf = 
        malloc(n_readers * sizeof(struct bam_reader_par));

    thread_params.reader_par = (void **)malloc(n_readers * sizeof(void *));

    unsigned r;
    unsigned long n_total_loci;
    if (query_range_file)
    {
        thread_params.ranges = 
            parse_query_ranges(query_range_file, &thread_params.n_ranges,
                               &n_total_loci);
        if (! thread_params.n_ranges)
        {
            fprintf(stderr, "Error: there are no ranges to process in query range file %s\n",
                    query_range_file);
            exit(1);
        }
    }
    else
    {
        /* simply set the 'query' to the largest span possible */
        thread_params.n_ranges = 1;
        thread_params.ranges = malloc(sizeof(struct pair_ordering_range));
        thread_params.ranges[0] = 
            (struct pair_ordering_range){ { 0, 0 }, { SIZE_MAX, SIZE_MAX - 1 } };
        n_total_loci = ULONG_MAX;
    }

    unsigned s;
    for (r = 0; r != n_readers; ++r)
    {
        thread_params.reader_buf[r] = (struct bam_reader_par){
            malloc(sample_attrs.n * sizeof(struct bam_stats)),
            sample_attrs.n,
            
            thread_params.ranges, 
            thread_params.ranges + thread_params.n_ranges
        };
        for (s = 0; s != sample_attrs.n; ++s)
        {
            thread_params.reader_buf[r].s[s].idx = NULL/* parse BAM index */;
            thread_params.reader_buf[r].s[s].bgzf = NULL/* open bam file */;
        }
        thread_params.reader_par[r] = &thread_params.reader_buf[r];
    }

    if (query_range_file)
        cs_init_by_range(n_total_loci, sample_attrs.n);

    else
    {
        cs_init_whole_file(sample_attrs.n);
        /* for (s = 0; s != sample_attrs.n; ++s) */
        /*     cs_set_total_bytes(s, ix[s].root->end_offset - ix[s].root->start_offset); */
    }

#define MAX_BYTES_SMALL_CHUNK 1e9
#define SMALL_CHUNK 5e6
#define DEFAULT_BYTES_PER_LOCUS 100

    cs_set_defaults(MAX_BYTES_SMALL_CHUNK,
                    SMALL_CHUNK, 
                    DEFAULT_BYTES_PER_LOCUS);

    thread_params.offload_par = 
        (struct dist_worker_offload_par){ dist_fh, comp_fh, indel_fh };

    /* To avoid a stall, n_extra / n_threads should be greater than
       Max(work chunk time) / Avg(work chunk time). */
    size_t n_extra = n_threads * NUM_EXTRA_BUFS_PER_THREAD;
    size_t n_output_files = 
        (dist_fh ? 1 : 0) + (comp_fh ? 1 : 0) + (indel_fh ? 1 : 0);

    thread_queue_reader_t reader = { bam_reader, bam_scanner };

    struct thread_queue *tqueue =
        thread_queue_init(reader, thread_params.reader_par,
                          dist_worker,
                          dist_offload, &thread_params.offload_par,
                          dist_on_create,
                          dist_on_exit,
                          n_threads, n_extra, n_readers, sample_attrs.n,
                          n_output_files, max_input_mem);

    return tqueue;
}


void dist_worker_tq_free()
{
    unsigned r, s;
    for (r = 0; r != thread_params.n_readers; ++r)
    {
        for (s = 0; sample_attrs.n; ++s)
            bam_reader_par_free(&thread_params.reader_buf[r].s[s]);

        free(thread_params.reader_buf[r].s);
    }

    free(thread_params.reader_buf);
    free(thread_params.reader_par);
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
 
void print_basecomp_quantiles(COMP_QV quantile_values,
                              const double *means,
                              size_t n_quantiles,
                              const char *label_string,
                              struct pileup_locus_info *pli,
                              struct pileup_data *pdat,
                              struct managed_buf *mb)
{
    char line_label[2048];

    sprintf(line_label,
            "%s\t%s\t%i\t%c\t%u\t%u",
            label_string, 
            pli->refname,
            pli->pos,
            pli->refbase,
            pdat->used_read_depth, 
            pdat->read_depth);

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
    for (d = 0; d != NUM_NUCS; ++d)
    {
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
void compute_wsq_dist(const double *points1,
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
void compute_dist_wquantiles(double *square_dist_buf,
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
   one if its pair index is equal to PSEUDO_SAMPLE */
void print_distance_quantiles(const char *contig,
                              size_t position,
                              char ref_base,
                              size_t pair_index,
                              struct locus_data *ldat,
                              double *dist_quantile_values,
                              struct managed_buf *buf)
{
    int s[] = { sample_pairs.p[pair_index].s1,
                sample_pairs.p[pair_index].s2 };

    unsigned space = (2 * MAX_LABEL_LEN) + 3 + 100 + (10 * MAX_NUM_QUANTILES);

    ALLOC_GROW_TYPED(buf->buf, buf->size + space, buf->alloc);

    buf->size +=
        sprintf(buf->buf + buf->size, 
                "%s\t%s\t%s\t%Zu\t%c", 
                sample_attrs.atts[s[0]].label,
                s[1] == PSEUDO_SAMPLE ? "REF" : sample_attrs.atts[s[1]].label,
                contig,
                position,
                ref_base);
    
    unsigned q;
    for (q = 0; q != n_quantiles; ++q)
        buf->size += sprintf(buf->buf + buf->size, "\t%7.4f",
                             dist_quantile_values[q]);

    if (worker_options.do_print_pileup)
    {
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

    for (pi = 0; pi != sample_pairs.n; ++pi)
    {
        int sp[] = { sample_pairs.p[pi].s1, sample_pairs.p[pi].s2 };

        ld[0] = &tls_dw.ldat[sp[0]];
        ld[1] = sp[1] == PSEUDO_SAMPLE ? &tls_dw.pseudo_sample : &tls_dw.ldat[sp[1]];
        
        tls_dw.metrics.total++;
        ++tls_dw.pair_stats[pi].total;

        /* load this particular pair of distp into the bep */
        for (i = 0; i != 2; ++i) {
            tls_dw.bep.dist[i] = &ld[i]->distp;
            if (! ld[i]->init.base_ct) {
                ld[i]->base_ct = pileup_basecall_stats(sp[i]);
                ld[i]->init.base_ct = 1;
            }
        }

        diff_state = 
            cached_dirichlet_diff(ld[0]->base_ct.ct, ld[1]->base_ct.ct, &tls_dw.bep,
                                  &cacheable, &cache_was_set);

        tls_dw.pair_stats[pi].dist_count[diff_state]++;
        tls_dw.pair_stats[pi].cacheable += cacheable;
        tls_dw.pair_stats[pi].cache_was_set += cache_was_set;

        tls_dw.metrics.cacheable += cacheable;
        tls_dw.metrics.cache_was_set += cache_was_set;

        if (diff_state == CHANGED)
        {
            /* Finish sampling and do full distance marginal estimation */
            for (i = 0; i != 2; ++i)
            {
                if (! ld[i]->init.bqs_ct) {
                    ld[i]->init.bqs_ct = 1;
                    pileup_bqs_stats(sp[i], &ld[i]->bqs_ct, &ld[i]->n_bqs_ct);
                }

                struct distrib_points *dst = tls_dw.bep.dist[i];
                struct points_gen_par *pgp = dst->pgen.points_gen_par;
                pgp->observed = ld[i]->bqs_ct;
                pgp->n_observed = ld[i]->n_bqs_ct;
                
                /* Generate all points */
                unsigned perm[] = { 0, 1, 2, 3 };
                update_points_gen_params(dst, ld[i]->base_ct.ct, perm);
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

            if (test_quantile_value > min_dirichlet_dist)
            {
                ++tls_dw.pair_stats[pi].confirmed_changed;
                ld[0]->confirmed_changed = 1;
                ld[1]->confirmed_changed = 1;

                if (out_buf)
                {
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


void print_indel_distance_quantiles(const char *contig,
                                    size_t position,
                                    size_t pair_index,
                                    double *dist_quantile_values,
                                    indel_event *events,
                                    size_t n_events,
                                    struct managed_buf *mb)
{
    int s1 = sample_pairs.p[pair_index].s1,
        s2 = sample_pairs.p[pair_index].s2;

    unsigned space = (2 * MAX_LABEL_LEN) + 3 + 100 + (10 * MAX_NUM_QUANTILES);
    ALLOC_GROW_TYPED(mb->buf, mb->size + space, mb->alloc);

    mb->size += sprintf(mb->buf + mb->size, 
                        "%s\t%s\t%s\t%Zu", 
                        samples.atts[s1].label,
                        s2 == PSEUDO_SAMPLE ? "REF" : samples.atts[s2].label,
                        contig, position);
    
    for (size_t q = 0; q != n_quantiles; ++q)
        mb->size += sprintf(mb->buf + mb->size, "\t%.4f", dist_quantile_values[q]);

    unsigned indel_space = n_events * (10 + 10);
    ALLOC_GROW_TYPED(mb->buf, mb->size + indel_space, mb->alloc);

    // now print the indel event summary
    indel_event *eb = events, *ee = eb + n_events;
    while (eb != ee) {
        mb->size +=
            sprintf(mb->buf + mb->size, "%c%i", (eb == events ? '\t' : ','), eb->count1);
        eb++;
    }
    eb = events;
    while (eb != ee) {
        mb->size += sprintf(mb->buf + mb->size, "%c%i", (eb == events ? '\t' : ','), eb->count2);
        eb++;
    }
    eb = events;
    while (eb != ee) {
        unsigned grow = 2 + (eb->seq ? strlen(eb->seq) : 0);
        ALLOC_GROW_TYPED(mb->buf, mb->size + grow, mb->alloc);

        mb->size += sprintf(mb->buf + mb->size, "%c%c%s", 
                            (eb == events ? '\t' : ','), 
                            (eb->seq ? (eb->is_insertion ? '+' : '-') : '@'),
                            (eb->seq ? eb->seq : ""));
        eb++;
    }

    locus_sampling 
        *ls1 = &dw->lslist[s1],
        *ls2 = s2 == PSEUDO_SAMPLE ? &dw->pseudo_sample : &dw->lslist[s2];

    if (worker_options.do_print_pileup) {
        unsigned extra_space = 
            ls1->locus.bases_raw.size
            + ls1->locus.quality_codes.size
            + ls2->locus.bases_raw.size
            + ls2->locus.quality_codes.size
            + 50;

        ALLOC_GROW_TYPED(mb->buf, mb->size + extra_space, mb->alloc);
        mb->size += sprintf(mb->buf + mb->size,
                            "\t%Zu\t%s\t%s\t%Zu\t%s\t%s",
                            ls1->locus.read_depth,
                            ls1->locus.bases_raw.buf,
                            ls1->locus.quality_codes.buf,
                            ls2->locus.read_depth,
                            ls2->locus.bases_raw.buf,
                            ls2->locus.quality_codes.buf);
    }
    mb->size += sprintf(mb->buf + mb->size, "\n");
}


// print out all next distance quantiles for indels
void indel_distance_quantiles_aux(struct managed_buf *buf)
{
    size_t n_events;

    for (size_t pi = 0; pi != sample_pairs.n; ++pi)
    {
        int s[] = { sample_pairs.p[pi].s1, sample_pairs.p[pi].s2 };
        
        struct locus_data ld[] = {
            &tls_dw->ldat[s[0]], 
            s[1] == PSEUDO_SAMPLE ? &tls_dw->pseudo_sample : &tls_dw->ldat[s[1]]
        };

        unsigned i;
        for (i = 0; i != 2; ++i) {
            if (! ld[i]->init.indel_ct) {
                ld[i]->init.indel_ct = 1;
                pileup_indel_stats(s[i], &ld[i]->indel_ct, ld[i]->n_indel_ct);
            }
        }
        if (ld[0]->n_indel_ct == 0 && ld[1]->n_indel_ct == 0
            && min_dirichlet_dist > 0) continue;
        
        /* traverse the sorted indel counts in tandem */
        unsigned max_n_events = ld[0]->n_indel_ct + ld[1]->n_indel_ct + 2;
        struct indel_pair_count *indel_pairs = 
            malloc(max_n_events * sizeof(struct indel_pair_count));

        struct indel_count *ic1 = , *ic2
        unsigned i;
        while (
        indel_event *all_events = count_indel_types(ls1, ls2, &n_events);

        if (n_events >= 2) 
        {
            // need at least some indels, otherwise these loci don't differ
            // there will always be at least one event, the reads themselves.  (though the count may be zero)
            double *alpha1 = new double[n_events], *alpha2 = new double[n_events];
            for (size_t c = 0; c != n_events; ++c) 
            {
                alpha1[c] = all_events[c].count1 + 1;
                alpha2[c] = all_events[c].count2 + 1;
            }
        
            size_t bufsize = n_events * max_sample_points;
            double *points1 = new double[bufsize], *p1 = points1, *pe1 = points1 + bufsize;
            double *points2 = new double[bufsize], *p2 = points2;

            while (p1 != pe1) 
            {
                gsl_ran_dirichlet(dw->randgen, n_events, alpha1, p1);
                gsl_ran_dirichlet(dw->randgen, n_events, alpha2, p2);
                p1 += n_events;
                p2 += n_events;
            }

            compute_square_dist(points1, points2, max_sample_points, n_events,
                                dw->square_dist_buf);

            double test_quantile = 1.0 - posterior_confidence, test_quantile_value;
            
            // compute distance quantiles
            compute_marginal_quantiles(dw->square_dist_buf,
                                       max_sample_points,
                                       1, /* one dimensional */
                                       0, /* use the first dimension */
                                       &test_quantile,
                                       1, /* evaluate only one quantile */
                                       &test_quantile_value);

            test_quantile_value = sqrt(test_quantile_value);
            if (test_quantile_value >= min_dirichlet_dist)
            {
                compute_marginal_quantiles(dw->square_dist_buf,
                                           max_sample_points,
                                           1, /* one dimensional */
                                           0, /* use the first dimension */
                                           quantiles,
                                           n_quantiles,
                                           dw->dist_quantile_values);
                unsigned q;
                for (q = 0; q != n_quantiles; ++q)
                    dw->dist_quantile_values[q] = 
                        sqrt(dw->dist_quantile_values[q]) * ONE_OVER_SQRT2;

                print_indel_distance_quantiles(ls1->locus.reference, 
                                               ls1->locus.position, pi, 
                                               dw->dist_quantile_values, 
                                               all_events, n_events, buf);
            }

            delete[] points1;
            delete[] points2;
            delete[] alpha1;
            delete[] alpha2;
        }
        delete[] all_events;
    }
}


void comp_quantiles_aux(struct managed_buf *comp_buf)
{
    unsigned s;
    struct pileup_locus_info ploc;
    pileup_current_info(&ploc);

    COMP_QV comp_quantile_values;
    double comp_means[NUM_NUCS];

    for (s = 0; s != sample_attrs.n; ++s) {
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
            if (! ld->init.sample_data)
            {
                ld->init.sample_data = 1;
                pileup_current_data(s, &ld->sample_data);

            print_basecomp_quantiles(comp_quantile_values,
                                     comp_means,
                                     n_quantiles,
                                     sample_attrs.atts[s].label,
                                     &ploc,
                                     &ld->sample_data,
                                     comp_buf);
        }
    }
}



/* receives a certain number of in_bufs and a certain number of
   out_bufs.
   
   there is one struct locus_data for each input.  it's current
   field points to the current line being processed. 'gs' is a single
   index indicating the sample with the lowest 'current' among all of
   them.  it is this position that must be fully processed before any
   sample_attrs may advance.
*/
void
dist_worker(const struct managed_buf *in_bufs,
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
    for (s = 0; s != sample_attrs.n; ++s)
    {
        bam_inflate(&in_bufs[s], &bam);
        tally_pileup_stats(bam, s);
    }

    /* */
    for (s = 0; s != sample_attrs.n; ++s)
        summarize_pileup_stats(s);

    free(bam.buf);

    /* zero out the pair_stats */
    unsigned pi;
    for (pi = 0; pi != sample_pairs.n; ++pi)
        memset(&tls_dw.pair_stats[pi], 0, sizeof(tls_dw.pair_stats[0]));

    
    unsigned more_loci = 1;
    while (more_loci)
    {
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

        more_loci = pileup_next_pos();
        for (s = 0; s != sample_attrs.n; ++s)
            reset_locus_data(tls_dw.ldat[s]);

    }   

    /* frees statistics that have already been used in one of the
       distance calculations. */
    pileup_clear_finished_stats();

    if (tls_dw.do_print_progress)
    {
        struct timespec now;
        clock_gettime(CLOCK_REALTIME, &now);
        unsigned elapsed = now.tv_sec - start_time.tv_sec;

        struct pileup_locus_info ploc;
        pileup_current_info(&ploc);

        time_t cal = time(NULL);
        char *ts = strdup(ctime(&cal));
        ts[strlen(ts)-1] = '\0';
        fprintf(stdout, 
                "%s (%02i:%02i:%02i elapsed): Finished processing %s %i\n", 
                ts,
                elapsed / 3600,
                (elapsed % 3600) / 60,
                elapsed % 60,
                ploc.refname,
                ploc.pos);
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
 
 
 void dist_offload(void *par, const struct managed_buf *bufs)
 {
     struct dist_worker_offload_par *ol = (struct dist_worker_offload_par *)par;
    unsigned i = 0;
    if (ol->dist_fh)
    {
        fwrite(bufs[i].buf, 1, bufs[i].size, ol->dist_fh), i++;
        fflush(ol->dist_fh);
    }

    if (ol->comp_fh)
    {
        fwrite(bufs[i].buf, 1, bufs[i].size, ol->comp_fh), i++;
        fflush(ol->comp_fh);
    }

    if (ol->indel_fh)
    {
        fwrite(bufs[i].buf, 1, bufs[i].size, ol->indel_fh), i++;
        fflush(ol->indel_fh);
    }
}
 
 
 
#if 0
/* update 'ls' fields to be consistent with sd->current */
void update_pileup_locus(const struct nucleotide_stats *stats,
                         struct locus_data *ls)
{
    ls->locus.load_line(ls->current);
    ls->locus_ord = init_locus(ls->current);
    ls->locus.parse();
    nucleotide_stats_pack(stats, &ls->locus.counts);
    ls->distp.points.size = 0;
    ls->distp.weights.size = 0;
}


void init_pseudo_locus(struct locus_data *ls)
{
    ls->is_next = 1;
    ls->locus.read_depth = PSEUDO_DEPTH;
    ls->locus.read_depth_match = PSEUDO_DEPTH;
    ls->locus.read_depth_high_qual = PSEUDO_DEPTH;

    ls->distp.points.size = 0;
    ls->distp.weights.size = 0;

    ((struct points_gen_par *)ls->distp.pgen.points_gen_par)->post_counts = 
        &ls->locus.counts;

    ls->locus.bases_raw.size = 4;
    ALLOC_GROW_TYPED(ls->locus.bases_raw.buf, 
                     ls->locus.bases_raw.size,
                     ls->locus.bases_raw.alloc);
    strcpy(ls->locus.bases_raw.buf, "REF");

    ls->locus.quality_codes.size = 4;
    ALLOC_GROW_TYPED(ls->locus.quality_codes.buf, 
                     ls->locus.quality_codes.size,
                     ls->locus.quality_codes.alloc);
    strcpy(ls->locus.quality_codes.buf, "REF");

}


/* update ls->locus.counts with ultra-high depth of perfect 'nuc'
   basecalls. */
void update_pseudo_locus(char nuc, struct locus_data *ls)
{
    unsigned inuc = (unsigned)base_to_index(nuc);
    if (inuc == 4)
        ls->locus.counts.num_data = 0;
    else
    {
        ls->locus.counts.num_data = 1;
        ls->locus.counts.stats_index[0] = 
            encode_nucleotide(nuc, NUC_HIGHEST_QUALITY - 1, 0);
        unsigned b;
        for (b = 0; b != NUM_NUCS; ++b)
        {
            ls->locus.counts.stats[0].cpd[b] = b == inuc ? 1 : 0;
            ls->locus.base_counts_high_qual[b] = b == inuc ? PSEUDO_DEPTH : 0;
        }

        ls->locus.counts.stats[0].ct = PSEUDO_DEPTH;
    }
}
#endif
