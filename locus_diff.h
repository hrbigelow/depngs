#ifndef LOCUS_DIFF_H
#define LOCUS_DIFF_H

#include <gsl/gsl_rng.h>

#include "defs.h"
#include "cache.h"
#include "thread_queue.h"
#include "dirichlet_diff_cache.h"
#include "dirichlet_points_gen.h"
#include "batch_pileup.h"
#include "bam_reader.h"

/* program-wide initialization of static variables */
void
locus_diff_init(double _post_confidence, 
                double _beta_confidence,
                double _min_dirichlet_dist,
                unsigned _max_sample_points,
                double _indel_prior_alpha,
                unsigned _max_dir_cache_items,
                unsigned _max_bounds_cache_items,
                unsigned n_threads,
                double prior_alpha,
                const char *samples_file,
                const char *sample_pairs_file,
                const char *fasta_file,
                struct bam_filter_params bam_filter,
                const char *quantiles_string,
                unsigned do_dist,
                unsigned do_comp,
                unsigned do_indel,
                unsigned do_print_pileup);

void locus_diff_free();


/* program-wide initialization of static variables, specific to the
   thread-queue. (also calls thread_queue_init) */
struct thread_queue *
locus_diff_tq_init(const char *query_range_file,
                   const char *fasta_file,
                   unsigned n_threads,
                   unsigned n_readers,
                   unsigned long max_input_mem,
                   FILE *dist_fh,
                   FILE *comp_fh,
                   FILE *indel_fh);

void locus_diff_tq_free();


/* there will be one of these instantiated for each thread.  Each of
   these holds the parameters needed by the thread that can be shared
   across different samples.  Since each thread computes a chunk of
   input across all samples, there are sample-specific parameters as
   well, held in sample_attributes. */
struct locus_diff_input
{
    struct binomial_est_params bep;
    size_t thread_num;

    double dist_quantile_values[MAX_NUM_QUANTILES];
    double comp_quantile_values[MAX_NUM_QUANTILES];

    struct {
        unsigned total, cacheable, cache_was_set;
    } metrics;

    struct locus_data pseudo_sample, *ldat;
    struct pair_dist_stats *pair_stats;
    gsl_rng *randgen;
    double *square_dist_buf; /* holds squares of distances for distance calculation */
    double *weights_buf; /* holds weights from those square distances
                            (product of weights on individual
                            points) */
    
    unsigned do_print_progress;
};


struct locus_diff_offload_par {
    FILE *dist_fh, *comp_fh, *indel_fh;
};

/* conforms to thread_queue_worker_t */
void
locus_diff_worker(const struct managed_buf *in_bufs,
                  unsigned more_input,
                  void *vsi, /* cast to struct bam_scanner_info */
                  struct managed_buf *out_bufs);

/* conforms to thread_queue_offload_t */
void locus_diff_offload(void *par, const struct managed_buf *bufs);

#define PSEUDO_DEPTH 1e6



void print_cache2_histo();

#endif // LOCUS_DIFF_H
