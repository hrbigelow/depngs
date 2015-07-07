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
#include "locus.h"
#include "dirichlet_points_gen.h"
#include "dirichlet_diff_cache.h"
#include "khash.h"
#include "range_line_reader.h"
#include "geometry.h"
#include "chunk_strategy.h"
#include "nucleotide_stats.h"
#include "common_tools.h"
#include "fastq_tools.h"
#include "dist_aux.h"

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
    struct range_line_reader_par *reader_buf;
    void **reader_par;
    unsigned n_readers;
    struct pair_ordering_range *ranges;
    unsigned n_ranges;
    unsigned n_threads;
    struct dist_worker_offload_par offload_par;
} thread_params;


__thread struct dist_worker_input tls_dw;


static double quantiles[MAX_NUM_QUANTILES];
static unsigned n_quantiles;




/*
void alloc_pileup_locus(struct locus_sampling *ls)
{
    ls->locus = PileupSummary();
    ls->is_next = 0;
    ls->confirmed_changed = 0;
    alloc_distrib_points(&ls->distp);
    ls->locus_ord = (struct pair_ordering){ 0, 0 };
}


void free_pileup_locus(struct locus_sampling *ls)
{
    free_distrib_points(&ls->distp);
}
*/


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
    /* alloc_pileup_locus(&tls_dw.pseudo_sample); */
    /* tls_dw.lslist = new struct locus_sampling[samples.n]; */
    /* for (s = 0; s != samples.n; ++s) */
    /*     alloc_pileup_locus(&tls_dw.lslist[s]); */

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
    /* free_pileup_locus(&tls_dw.pseudo_sample); */
    /* for (s = 0; s != samples.n; ++s) */
    /*     free_pileup_locus(&tls_dw.lslist[s]); */

    /* delete[] tls_dw.lslist; */
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

    size_t scan_thresh_size = 1e5;
    file_bsearch_init(init_locus, scan_thresh_size);

    thread_params.reader_buf = 
        malloc(n_readers * sizeof(struct range_line_reader_par));

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
        thread_params.reader_buf[r] = (struct range_line_reader_par){
            malloc(sample_attrs.n * sizeof(struct file_bsearch_index)),
            sample_attrs.n,
            thread_params.ranges, 
            thread_params.ranges + thread_params.n_ranges,
            NULL,
            0,
            init_locus,
            1
        };
        for (s = 0; s != sample_attrs.n; ++s)
            thread_params.reader_buf[r].ix[s] =
                file_bsearch_make_index(sample_attrs.atts[s].file);

        thread_params.reader_par[r] = &thread_params.reader_buf[r];
    }

    if (query_range_file)
        cs_init_by_range(n_total_loci, sample_attrs.n);
    else
    {
        cs_init_whole_file(sample_attrs.n);
        struct file_bsearch_index *ix = thread_params.reader_buf[0].ix;
        for (s = 0; s != sample_attrs.n; ++s)
            cs_set_total_bytes(s, ix[s].root->end_offset - ix[s].root->start_offset);
    }
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

    thread_queue_reader_t reader = { rl_reader, rl_scanner };

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
            file_bsearch_index_free(thread_params.reader_buf[r].ix[s]);

        free(thread_params.reader_buf[r].ix);
    }

    free(thread_params.reader_buf);
    free(thread_params.reader_par);
    free(thread_params.ranges);

}

typedef double COMP_QV[NUM_NUCS][MAX_NUM_QUANTILES];

void print_basecomp_quantiles(COMP_QV quantile_values,
                              const double *means,
                              size_t n_quantiles,
                              const char *label_string,
                              struct locus_sampling *ls, 
                              struct managed_buf *mb)
{
    char line_label[2048];

    sprintf(line_label,
            "%s\t%s\t%i\t%c\t%Zu\t%Zu",
            label_string, 
            ls->locus.reference, 
            ls->locus.position, 
            ls->locus.reference_base, 
            ls->locus.read_depth, 
            ls->locus.read_depth_high_qual
            );

    std::multimap<double, size_t, std::greater<double> > dim_to_mean;
    unsigned d, q;
    for (d = 0; d != NUM_NUCS; ++d)
        dim_to_mean.insert(std::make_pair(means[d], d));

    /* calculate mean rank order */
    size_t mean_rank_order[NUM_NUCS];
    std::multimap<double, size_t, std::greater<double> >::iterator dtm;
    d = 0;
    for (dtm = dim_to_mean.begin(); dtm != dim_to_mean.end(); ++dtm)
        mean_rank_order[(*dtm).second] = d++;

    unsigned grow = NUM_NUCS * (sizeof(line_label) + 30 + (11 * MAX_NUM_QUANTILES));
    ALLOC_GROW_TYPED(mb->buf, mb->size + grow, mb->alloc);

    static const char dimension_labels[] = "ACGT";
    for (d = 0; d != NUM_NUCS; ++d)
    {
        mb->size += sprintf(mb->buf + mb->size,
                            "%s\t%c\t%Zu\t%10.8f", line_label, 
                            dimension_labels[d], mean_rank_order[d], means[d]);
        
        for (q = 0; q != n_quantiles; ++q)
            mb->size += sprintf(mb->buf + mb->size, "\t%10.8f", quantile_values[d][q]);

        mb->size += sprintf(mb->buf + mb->size, "\n");
    }
}




/* populates square_dist_buf with squares of euclidean distances
   between points1 and points2 in the barycentric space (R4,
   normalized positive components).  populates weights_buf with
   product of weights1 and weights2. */
void compute_wsq_dist(const POINT *points1,
                      const double *weights1,
                      const POINT *points2,
                      const double *weights2,
                      size_t n_points,
                      double *square_dist_buf,
                      double *weights_buf)
{
    compute_square_dist((const double *)points1, (const double *)points2, n_points, 4, square_dist_buf);
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
                              double *dist_quantile_values,
                              struct managed_buf *buf)
{
    int s1 = sample_pairs.p[pair_index].s1,
        s2 = sample_pairs.p[pair_index].s2;

    unsigned space = (2 * MAX_LABEL_LEN) + 3 + 100 + (10 * MAX_NUM_QUANTILES);

    ALLOC_GROW_TYPED(buf->buf, buf->size + space, buf->alloc);

    buf->size += sprintf(buf->buf + buf->size, 
                         "%s\t%s\t%s\t%Zu\t%c", 
                         sample_attrs.atts[s1].label,
                         s2 == PSEUDO_SAMPLE ? "REF" : sample_attrs.atts[s2].label,
                         contig,
                         position,
                         ref_base);
    
    for (size_t q = 0; q != n_quantiles; ++q)
        buf->size += sprintf(buf->buf + buf->size, "\t%7.4f", dist_quantile_values[q]);

    locus_sampling *ls1 = &tls_dw.lslist[s1], 
        *ls2 = s2 == PSEUDO_SAMPLE ? &tls_dw.pseudo_sample : &tls_dw.lslist[s2];

    if (worker_options.do_print_pileup)
    {
        unsigned extra_space = 
            ls1->locus.bases_raw.size
            + ls1->locus.quality_codes.size
            + ls2->locus.bases_raw.size
            + ls2->locus.quality_codes.size
            + 50;

        ALLOC_GROW_TYPED(buf->buf, buf->size + extra_space, buf->alloc);
        
        buf->size += sprintf(buf->buf + buf->size, 
                             "\t%Zu\t%s\t%s\t%Zu\t%s\t%s",
                             ls1->locus.read_depth,
                             ls1->locus.bases_raw.buf,
                             ls1->locus.quality_codes.buf,
                             ls2->locus.read_depth,
                             ls2->locus.bases_raw.buf,
                             ls2->locus.quality_codes.buf);
    }
    buf->size += sprintf(buf->buf + buf->size, "\n");
}

/* update 'ls' fields to be consistent with sd->current */
void update_pileup_locus(const struct nucleotide_stats *stats,
                         locus_sampling *ls)
{
    ls->locus.load_line(ls->current);
    ls->locus_ord = init_locus(ls->current);
    ls->locus.parse();
    nucleotide_stats_pack(stats, &ls->locus.counts);
    ls->distp.points.size = 0;
    ls->distp.weights.size = 0;
}


void init_pseudo_locus(struct locus_sampling *ls)
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
void update_pseudo_locus(char nuc, struct locus_sampling *ls)
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

/* advance ls->current, then initialize if there is another locus */
void refresh_locus(const struct nucleotide_stats *stats,
                   locus_sampling *ls)
{
    assert(ls->current != ls->end);
    if ((ls->current = strchr(ls->current, '\n') + 1) == ls->end) 
        ls->is_next = 0;

    else
    {
        update_pileup_locus(stats, ls);
        ls->confirmed_changed = 0;
    }
}



#define ONE_OVER_SQRT2 0.70710678118654752440

/* for each sample pair, calculate whether the loci differ above the
   given level and confidence. if a pair differs, set the
   confirmed_changed flag for each sample in the pair. if out_buf is
   not NULL, print out distance quantiles.  also generates sample
   points for each sample as needed, both for the preliminary test and
   more points for the final test */
void next_distance_quantiles_aux(size_t gs,
                                 struct managed_buf *out_buf)
{
    locus_sampling *ls1, *ls2;
    size_t pi, i;
    unsigned counts[2][NUM_NUCS];
    enum fuzzy_state diff_state = AMBIGUOUS;

    unsigned cacheable, cache_was_set;

    for (pi = 0; pi != sample_pairs.n; ++pi)
    {
        int s1 = sample_pairs.p[pi].s1,
            s2 = sample_pairs.p[pi].s2;

        ls1 = &tls_dw.lslist[s1];
        ls2 = s2 == PSEUDO_SAMPLE ? &tls_dw.pseudo_sample : &tls_dw.lslist[s2];
        
        tls_dw.bep.dist[0] = &ls1->distp;
        tls_dw.bep.dist[1] = &ls2->distp;

        for (i = 0; i != NUM_NUCS; ++i)
        {
            counts[0][i] = ls1->locus.base_counts_high_qual[i];
            counts[1][i] = ls2->locus.base_counts_high_qual[i];
        }

        tls_dw.metrics.total++;

        ++tls_dw.pair_stats[pi].total;

        if (! (ls1->is_next && ls2->is_next))
            diff_state = AMBIGUOUS, cacheable = 1, cache_was_set = 1;

        else
            diff_state = 
                cached_dirichlet_diff(counts[0], counts[1], &tls_dw.bep,
                                      &cacheable, &cache_was_set);

        ++tls_dw.pair_stats[pi].dist_count[diff_state];
        tls_dw.pair_stats[pi].cacheable += cacheable;
        tls_dw.pair_stats[pi].cache_was_set += cache_was_set;

        tls_dw.metrics.cacheable += cacheable;
        tls_dw.metrics.cache_was_set += cache_was_set;

        if (diff_state == CHANGED)
        {
            /* Finish sampling and do full distance marginal estimation */
            for (i = 0; i != 2; ++i)
            {
                unsigned perm[] = { 0, 1, 2, 3 };
                struct distrib_points *dst = tls_dw.bep.dist[i];
                /* Generate all points */
                update_points_gen_params(dst, counts[i], perm);
                POINT *p,
                    *pb = dst->points.buf,
                    *pe = dst->points.buf + max_sample_points;
                for (p = pb; p != pe; p += GEN_POINTS_BATCH)
                    dst->pgen.gen_point(dst->pgen.points_gen_par, p);

                dst->points.size = max_sample_points;
                
                /* Generate all weights */
                double *w,
                    *wb = dst->weights.buf,
                    *we = dst->weights.buf + max_sample_points;
                for (w = wb, p = pb; w != we; 
                     w += GEN_POINTS_BATCH, p += GEN_POINTS_BATCH)
                    dst->pgen.weight(p, dst->pgen.points_gen_par, w);

                dst->weights.size = max_sample_points;
            }
            /* Compute weighted square distances (max value of 2).
               e.g.
               minimum simplex distance on [0,1.00] scale is 0.25.
               minimum real    distance on [0,1.41] scale is 0.35355
               minimum sq_real distance on [0,2.00] scale is 0.12499
            */

            compute_wsq_dist(tls_dw.bep.dist[0]->points.buf, 
                             tls_dw.bep.dist[0]->weights.buf,
                             tls_dw.bep.dist[1]->points.buf, 
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
                ls1->confirmed_changed = 1;
                ls2->confirmed_changed = 1;

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
                    
                    print_distance_quantiles(tls_dw.lslist[gs].locus.reference, 
                                             tls_dw.lslist[gs].locus.position, 
                                             tls_dw.lslist[gs].locus.reference_base,
                                             pi, tls_dw.dist_quantile_values, 
                                             out_buf);
                }
            }
        }
        
    }
}




// find gs, and init the is_next field for all sample_attrs
// this is run
size_t init_global_and_next(locus_sampling *lslist)
{
    size_t gs = 0;
    for (size_t s = 0; s != sample_attrs.n; ++s)
    {
        if (lslist[s].current != lslist[s].end
            && (lslist[gs].current == lslist[gs].end
                || cmp_pair_ordering(&lslist[s].locus_ord, 
                                     &lslist[gs].locus_ord) < 0)
            )
            gs = s;
    }
    bool global_at_end = (lslist[gs].current == lslist[gs].end);
    unsigned s;
    for (s = 0; s != sample_attrs.n; ++s)
        lslist[s].is_next = 
            lslist[s].current != lslist[s].end
            && (! global_at_end)
            && cmp_pair_ordering(&lslist[s].locus_ord,
                                 &lslist[gs].locus_ord) == 0;

    return gs;
}


/* receives a certain number of in_bufs and a certain number of
   out_bufs.
   
   there is one struct locus_sampling for each input.  it's current
   field points to the current line being processed. 'gs' is a single
   index indicating the sample with the lowest 'current' among all of
   them.  it is this position that must be fully processed before any
   sample_attrs may advance.

   any sample missing the locus defined by sample[gs].current has
   null_sd substituted for it.
 */
void dist_worker(const struct managed_buf *in_bufs,
                 struct managed_buf *out_bufs)
{
    struct timespec worker_start_time;
    clock_gettime(CLOCK_REALTIME, &worker_start_time);
    tls_dw.metrics = { 0, 0, 0 };

    unsigned i = 0;
    struct managed_buf 
        *dist_buf = worker_options.do_dist ? &out_bufs[i++] : NULL,
        *comp_buf = worker_options.do_comp ? &out_bufs[i++] : NULL,
        *indel_buf = worker_options.do_indel ? &out_bufs[i++] : NULL;

    size_t gs = 0, s;
    struct locus_sampling *ls;
    for (s = 0; s != sample_attrs.n; ++s)
    {    
        ls = &tls_dw.lslist[s];
        ls->current = in_bufs[s].buf;
        ls->end = in_bufs[s].buf + in_bufs[s].size;
        ((struct points_gen_par *)ls->distp.pgen.points_gen_par)->post_counts = 
            &ls->locus.counts;
    }

    // before main loop, initialize loci and sample_attrs
    for (s = 0; s != sample_attrs.n; ++s)
    {
        if (tls_dw.lslist[s].current == tls_dw.lslist[s].end) continue;
        update_pileup_locus(&sample_attrs.atts[s].nuc_stats, &tls_dw.lslist[s]);
    }

    init_pseudo_locus(&tls_dw.pseudo_sample);

    gs = init_global_and_next(tls_dw.lslist);

    update_pseudo_locus(tls_dw.lslist[gs].locus.reference_base, &tls_dw.pseudo_sample);

    /* zero out the pair_stats */
    unsigned pi;
    for (pi = 0; pi != sample_pairs.n; ++pi)
        memset(&tls_dw.pair_stats[pi], 0, sizeof(tls_dw.pair_stats[0]));

    // main loop for computing pairwise distances
    // though slightly inefficient, just construct null_points here
    // !!! here, use worker[0] as a proxy.  this is a design flaw
    // owing to the fact that there are multiple workers used, all with the same
    // basic settings
    locus_sampling null_sd;
    alloc_pileup_locus(&null_sd);

    char null_pileup[50];
    strcpy(null_pileup, tls_dw.lslist[0].locus.reference);
    strcat(null_pileup, "\t1\tA\t0\t\t\n");
    null_sd.current = null_pileup;
    null_sd.end = null_pileup + strlen(null_pileup);
        
    update_pileup_locus(&sample_attrs.atts[0].nuc_stats, &null_sd);

    COMP_QV comp_quantile_values;
    double comp_means[NUM_NUCS];

    while (tls_dw.lslist[gs].current != tls_dw.lslist[gs].end)
    {
        /* this will strain the global mutex, but only during the hash
           loading phase. */
        if (! tls_dw.bep.points_hash_frozen)
            tls_dw.bep.points_hash_frozen = freeze_points_hash();

        if (! tls_dw.bep.bounds_hash_frozen)
            tls_dw.bep.bounds_hash_frozen = freeze_bounds_hash();

        if (dist_buf || comp_buf)
            next_distance_quantiles_aux(gs, dist_buf);

        if (indel_buf)
            next_indel_distance_quantiles_aux(gs, indel_buf);

        if (comp_buf)
        {
            for (s = 0; s != sample_attrs.n; ++s)
            {
                struct locus_sampling *sam = &tls_dw.lslist[s];
                if (sam->is_next && sam->confirmed_changed)
                {
                    unsigned d;
                    for (d = 0; d != NUM_NUCS; ++d)
                    {
                        compute_marginal_wquantiles((double *)sam->distp.points.buf,
                                                    sam->distp.weights.buf,
                                                    max_sample_points,
                                                    NUM_NUCS,
                                                    d,
                                                    quantiles,
                                                    n_quantiles,
                                                    comp_quantile_values[d]);
                        comp_means[d] = 
                            compute_marginal_mean((double *)sam->distp.points.buf,
                                                  sam->distp.weights.buf,
                                                  max_sample_points,
                                                  NUM_NUCS,
                                                  d);
                    }
                    
                    print_basecomp_quantiles(comp_quantile_values,
                                             comp_means,
                                             n_quantiles,
                                             sample_attrs.atts[s].label,
                                             &tls_dw.lslist[s],
                                             comp_buf);
                }
            }
        }
        for (s = 0; s != sample_attrs.n; ++s)
            if (tls_dw.lslist[s].is_next) 
                refresh_locus(&sample_attrs.atts[s].nuc_stats, &tls_dw.lslist[s]);

        gs = init_global_and_next(tls_dw.lslist);

        /* refresh the pseudo sample */
        update_pseudo_locus(tls_dw.lslist[gs].locus.reference_base, &tls_dw.pseudo_sample);
    }   

    if (tls_dw.do_print_progress)
    {
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
                tls_dw.lslist[gs].locus.reference,
                tls_dw.lslist[gs].locus.position);
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

    free_pileup_locus(&null_sd);
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


