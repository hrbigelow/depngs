#ifndef _BINOMIAL_EST_H
#define _BINOMIAL_EST_H

#include "defs.h"
#include <stdlib.h>

enum fuzzy_state {
    CHANGED,                /* the two loci differ */
    AMBIGUOUS_OR_CHANGED,   /* max_points taken; UNCHANGED category eliminated  */
    AMBIGUOUS,              /* posterior distributions too diffuse to call */
    AMBIGUOUS_OR_UNCHANGED, /* max_points taken; CHANGED category eliminated */
    UNCHANGED               /* the two loci do not differ */
};


struct binomial_est_state {
    enum fuzzy_state state;
    double beta_qval_lo, beta_qval_hi;
};


typedef double POINT[NUM_NUCS];

struct points_buf {
    POINT *buf;
    size_t size, alloc;
};

struct weights_buf {
    double *buf;
    size_t size, alloc;
};

struct points_gen
{
    void *point_par;
    void (*gen_point)(const void *par, POINT *points);

    void *weight_par;
    void (*weight)(POINT *points, const void *par,
                   double *weights);
};


void init_beta(double beta_conf);


/* Sample pairs of points from dist_pair up to max_points, classifying
   each pair as 'success' if distance is less than min_dist, 'failure'
   otherwise.  From the set of successes and failures, use the Beta
   distribution to estimate the true binomial probability.  Use dist1
   and dist2 to generate more points as needed. */
struct binomial_est_state
binomial_quantile_est(unsigned max_points, float min_dist,
                      float post_conf, float beta_conf,
                      struct points_gen pgen1,
                      struct points_buf *points1,
                      struct points_gen pgen2,
                      struct points_buf *points2,
                      size_t batch_size);

#endif /* _BINOMIAL_EST_H */
