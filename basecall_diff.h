#ifndef DIST_WORKER_H
#define DIST_WORKER_H

#include <gsl/gsl_rng.h>

#include "defs.h"
#include "cache.h"
// #include "pileup_tools.h"

/* #include "binomial_est.h" */
/* #include "ordering.h" */
#include "thread_queue.h"
#include "dirichlet_diff_cache.h"

/* instantiate one of these for each thread and each sample.  holds
   the information for the locus currently being processed in this
   thread and sample. */
/*
struct locus_sampling
{
    PileupSummary locus;
    unsigned char is_next;
    unsigned char confirmed_changed;
    struct distrib_points distp;
    char *current, *end;
    pair_ordering locus_ord;
};
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
                      unsigned do_print_pileup);

void dist_worker_free();

struct thread_queue *
dist_worker_tq_init(const char *query_range_file,
                    unsigned n_threads,
                    unsigned n_readers,
                    unsigned long max_input_mem,
                    FILE *dist_fh,
                    FILE *comp_fh,
                    FILE *indel_fh);

void dist_worker_tq_free();



/* there will be one of these instantiated for each thread.  Each of
   these holds the parameters needed by the thread that can be shared
   across different samples.  Since each thread computes a chunk of
   input across all samples, there are sample-specific parameters as
   well, held in sample_attributes. */
struct dist_worker_input
{
    struct binomial_est_params bep;
    size_t thread_num;

    double dist_quantile_values[MAX_NUM_QUANTILES];
    double comp_quantile_values[MAX_NUM_QUANTILES];

    struct {
        unsigned total, cacheable, cache_was_set;
    } metrics;

    /* struct locus_sampling pseudo_sample, *lslist; */
    struct pair_dist_stats *pair_stats;
    gsl_rng *randgen;
    double *square_dist_buf; /* holds squares of distances for distance calculation */
    double *weights_buf; /* holds weights from those square distances
                            (product of weights on individual
                            points) */

    unsigned do_print_progress;
};


struct dist_worker_offload_par {
    FILE *dist_fh, *comp_fh, *indel_fh;
};

/* conforms to thread_queue_worker_t */
void dist_worker(void *par,
                 const struct managed_buf *in_bufs,
                 struct managed_buf *out_bufs);

/* conforms to thread_queue_offload_t */
void dist_offload(void *par, const struct managed_buf *bufs);

#define PSEUDO_DEPTH 100000

void print_cache2_histo();

#endif // DIST_WORKER_H
