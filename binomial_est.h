#ifndef _BINOMIAL_EST_H
#define _BINOMIAL_EST_H

#include "defs.h"
#include "dirichlet_points_gen.h"
#include <stdlib.h>

struct binomial_est_params {
    unsigned max_sample_points;
    double post_confidence;
    double beta_confidence;
    double min_dirichlet_dist;
    unsigned batch_size;
};

enum fuzzy_state {
    CHANGED,                /* the two loci differ */
    AMBIGUOUS_OR_CHANGED,   /* max_points taken; UNCHANGED category eliminated  */
    AMBIGUOUS,              /* posterior distributions too diffuse to call */
    AMBIGUOUS_OR_UNCHANGED, /* max_points taken; CHANGED category eliminated */
    UNCHANGED               /* the two loci do not differ */
};


extern const char *fuzzy_state_strings[];

struct binomial_est_state {
    enum fuzzy_state state;
    double beta_qval_lo, beta_qval_hi;
};


struct binomial_est_bounds {
    /* enum init_phase state; */
    int32_t ambiguous[2];
    int32_t unchanged[2];
};





/* initializes a private global beta_cache */
void binomial_est_init(struct binomial_est_params be_par,
                       unsigned num_beta_precalc,
                       size_t n_threads);


/* frees the beta_cache */
void binomial_est_free();


/* Sample pairs of points from dist_pair up to max_points, classifying
   each pair as 'success' if distance is less than min_dist, 'failure'
   otherwise.  From the set of successes and failures, use the Beta
   distribution to estimate the true binomial probability.  Use dist1
   and dist2 to generate more points as needed. */
struct binomial_est_state
binomial_quantile_est(struct points_gen pgen1,
                      struct points_buf *points1,
                      struct points_gen pgen2,
                      struct points_buf *points2);

#endif /* _BINOMIAL_EST_H */
