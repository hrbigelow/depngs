#ifndef LOCUS_DIFF_H
#define LOCUS_DIFF_H

#include <gsl/gsl_rng.h>

#include "defs.h"
#include "bam_reader.h"
#include "dir_cache.h"
#include "dirichlet_diff_cache.h"
#include "thread_queue.h"

struct locus_diff_params {
    unsigned do_print_pileup;
    unsigned do_dist, do_comp, do_indel;
    unsigned max_sample_points;
    double post_confidence;
    double min_dirichlet_dist;
    double prior_alpha;
    double indel_prior_alpha;
    double quantiles[MAX_NUM_QUANTILES];
    unsigned n_quantiles;
};




/* program-wide initialization of static variables */
struct thread_queue *
locus_diff_init(const char *samples_file, const char *sample_pair_file, 
                const char *locus_range_file, const char *fasta_file,
                unsigned n_threads, unsigned n_max_reading, unsigned long max_input_mem,
                struct locus_diff_params ld_par,
                struct dirichlet_diff_params dd_par,
                struct binomial_est_params be_par,
                struct dir_cache_params dc_par,
                struct bam_filter_params bf_par,
                FILE *dist_fh, FILE *comp_fh, FILE *indel_fh);


void locus_diff_free();


/* program-wide initialization of static variables, specific to the
   thread-queue. (also calls thread_queue_init) */
struct thread_queue *
locus_diff_tq_init(const char *locus_range_file,
                   const char *fasta_file,
                   unsigned n_threads,
                   unsigned n_max_reading,
                   unsigned long max_input_mem,
                   struct dirichlet_diff_params dd_par,
                   struct binomial_est_params be_par,
                   struct dir_cache_params dc_par, 
                   struct bam_filter_params bf_par,
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
    struct bound_search_params bep;
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
