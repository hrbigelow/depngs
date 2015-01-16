#ifndef DIST_WORKER_H
#define DIST_WORKER_H

#include <map>
#include <vector>
#include <cstdlib>
#include <cstring>

#include <gsl/gsl_randist.h>

#include "locus_comp.h"
#include "defs.h"

struct sample_attributes;

/* there will be one of these instantiated for each thread.  Each of
   these holds the parameters needed by the thread that can be shared
   across different samples.  Since each thread computes a chunk of
   input across all samples, there are sample-specific parameters as
   well, held in sample_attributes. */
struct dist_worker_input
{
    struct sample_attributes *sample_atts;
    struct posterior_settings *pset;
    size_t thread_num;
    size_t n_samples;
    size_t n_sample_pairs;
    size_t n_sample_point_pairings; 

    double dist_quantile_values[MAX_NUM_QUANTILES];
    double comp_quantile_values[MAX_NUM_QUANTILES];

    float min_high_conf_dist; 
    size_t prelim_n_points;
    double prelim_quantile;

    gsl_rng *randgen;
    double *square_dist_buf; /* holds squares of distances for distance calculation */
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
