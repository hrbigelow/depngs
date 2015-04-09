#ifndef DIST_WORKER_H
#define DIST_WORKER_H

#include <map>
#include <vector>
#include <cstdlib>
#include <cstring>

#include <gsl/gsl_rng.h>

#include "defs.h"
#include "pileup_tools.h"

extern "C" {
#include "binomial_est.h"
#include "ordering.h"
#include "metropolis_sampling.h"
#include "thread_queue.h"
#include "dirichlet_diff_cache.h"
}

/* instantiate one of these for each thread and each sample.  holds
   the information for the locus currently being processed in this
   thread and sample. */
struct locus_sampling
{
    PileupSummary locus;
    bool is_next;
    unsigned char dist_printed;
    struct distrib_points distp;
    char *current, *end;
    pair_ordering locus_ord;
};

/* attributes intrinsic to one sample */
struct sample_attributes
{
    char label_string[100];
    struct nucleotide_stats nuc_stats;
    FILE *fh;
};

void init_sample_attributes(const char *jpd_file,
                            const char *sample_label,
                            const char *pileup_file,
                            struct sample_attributes *s);

struct pair_dist_stats {
    size_t dist_count[5]; /* number of enum fuzzy_state choices */
    size_t confirmed_changed;
};


/* there will be one of these instantiated for each thread.  Each of
   these holds the parameters needed by the thread that can be shared
   across different samples.  Since each thread computes a chunk of
   input across all samples, there are sample-specific parameters as
   well, held in sample_attributes. */
struct dist_worker_input
{
    struct sample_attributes *sample_atts;
    struct binomial_est_params bep;
    /* struct posterior_settings *pset; */
    size_t thread_num;
    size_t n_samples;
    size_t n_sample_pairs;
    size_t n_sample_point_pairings; 

    double dist_quantile_values[MAX_NUM_QUANTILES];
    double comp_quantile_values[MAX_NUM_QUANTILES];

    gsl_rng *randgen;
    double *square_dist_buf; /* holds squares of distances for distance calculation */
    double *weights_buf; /* holds weights from those square distances
                            (product of weights on individual
                            points) */

    struct pair_dist_stats *pair_stats;

    int print_pileup_fields; // if 0, do not print extra pileup fields

    int do_dist, do_comp, do_indel;

    /* defines the parsed set of sample pairs to compare */
    size_t *pair_sample1, *pair_sample2; 
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

#endif // DIST_WORKER_H
