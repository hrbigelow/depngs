#ifndef _DIRICHLET_DIFF_CACHE_H
#define _DIRICHLET_DIFF_CACHE_H

#include "binomial_est.h"
#include "metropolis_sampling.h"

struct distrib_points {
    struct points_gen pgen;
    struct points_buf points;
    struct weights_buf weights;
};


/* one instance per thread. */
struct binomial_est_params {
    enum fuzzy_state query_state;
    unsigned use_low_beta;
    double query_beta;
    struct posterior_settings *pset;
    struct distrib_points *dist[2];
    size_t batch_size;
};


void dirichlet_diff_init();
void dirichlet_diff_free();


void print_beb_bounds(struct binomial_est_params *bpar);

void print_bounds(struct binomial_est_params *bpar);


/* test two dirichlets based on their counts. Use a thread-safe
   caching mechanism to update the cache as necessary. */
enum fuzzy_state
cached_dirichlet_diff(unsigned *a_counts,
                      unsigned *b_counts,
                      struct binomial_est_params *bpar);


#endif /* _DIRICHLET_DIFF_CACHE_H */
