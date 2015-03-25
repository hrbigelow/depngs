#ifndef _METROPOLIS_SAMPLING_H
#define _METROPOLIS_SAMPLING_H

#include "nucleotide_stats.h"
#include "defs.h"

/* should be just one instance of this for an entire run.  holds
   settings for doing metropolis sampling, plus settings common to
   comp or dist estimation. */
struct posterior_settings
{
    double prior_alpha[NUM_NUCS];
    size_t max_sample_points; /* maximum number of sample points to take */
    double min_dist, post_confidence, beta_confidence;
    double min_quality_score;
    double dist_quantiles[MAX_NUM_QUANTILES];
    double comp_quantiles[MAX_NUM_QUANTILES];
    size_t n_dist_quantiles, n_comp_quantiles;
};

struct eval_counts {
    unsigned n_dirichlet;
    unsigned n_yep;
    unsigned n_tuning_iter;
    unsigned cumul_aoff;
};


double alphas_from_counts(const struct packed_counts *cts, 
                          const double *prior_alpha,
                          double *est_alpha);

size_t tune_proposal(const struct packed_counts *cts,
                     const struct posterior_settings *set,
                     double *proposal_alpha,
                     double *estimated_mean,
                     double *points_buf,
                     struct eval_counts *eval);

void
metropolis_sampling(unsigned short start_point, unsigned short n_points_wanted,
                    const struct packed_counts *cts,
                    const double *logu, /* logs of U[0, 1] values */
                    double *proposal_alpha,
                    const double *prior_alpha,
                    size_t nth, /* collect every nth point */
                    double *sample_points,
                    struct eval_counts *eval);

#endif /* METROPOLIS_SAMPLING_H */
