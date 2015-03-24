#ifndef _DIRICHLET_DIFF_CACHE_H
#define _DIRICHLET_DIFF_CACHE_H

#include "binomial_est.h"

/* Should have one of these per thread */
struct binomial_est_params {
    enum fuzzy_state query_state;
    unsigned max_points;
    float min_dist;
    float post_conf;
    float beta_conf;
    double prior_alpha[NUM_NUCS];
    struct points_gen pgen1;
    struct points_buf *points1;
    struct points_gen pgen2;
    struct points_buf *points2;
    size_t batch_size;
};


/* test two dirichlets based on their counts. Use a thread-safe
   caching mechanism to update the cache as necessary. */
enum fuzzy_state cached_dirichlet_diff(unsigned *a_counts,
                                       unsigned *b_counts,
                                       struct binomial_est_params *bpar);


#endif /* _DIRICHLET_DIFF_CACHE_H */
