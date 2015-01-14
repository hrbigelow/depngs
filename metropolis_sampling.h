#ifndef _METROPOLIS_SAMPLING_H
#define _METROPOLIS_SAMPLING_H

#include "nucleotide_stats.h"
#include "defs.h"

/* should be just one instance of this for an entire run */
struct posterior_settings
{
    double prior_alpha[4];
    size_t max_tuning_iterations;
    size_t tuning_n_points;
    size_t final_n_points;
    double autocor_max_offset;
    double autocor_max_value;

    size_t initial_autocor_offset;
    size_t target_autocor_offset;
    double *logu; /* log(U) for a set of values U sampled from Uniform(0, 1) */
    double min_quality_score;
    double dist_quantiles[MAX_NUM_QUANTILES];
    double comp_quantiles[MAX_NUM_QUANTILES];
    size_t n_dist_quantiles, n_comp_quantiles;
};


size_t tune_proposal(struct packed_counts *cts,
                     const struct posterior_settings *set,
                     double *proposal_alpha,
                     double *estimated_mean,
                     double *points_buf);

void metropolis_sampling(unsigned short start_point, unsigned short n_points_wanted,
                         const struct packed_counts *cts,
                         const double *logu, /* logs of U[0, 1] values */
                         double *proposal_alpha,
                         size_t nth, /* collect every nth point */
                         double *sample_points);

#endif /* METROPOLIS_SAMPLING_H */
