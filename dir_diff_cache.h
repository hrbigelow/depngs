#ifndef _DIR_DIFF_CACHE_H
#define _DIR_DIFF_CACHE_H

#include "binomial_est.h"
#include "dir_cache.h"

#include <stdio.h>

/* global configuration parameters for dirichlet_diff_cache */
struct dirichlet_diff_params {
    unsigned pseudo_depth;
    double prior_alpha;
} g_dd_par;


/* one instance per thread. */
struct bound_search_params {
    unsigned use_low_beta;
    double query_beta;
    /* these two are set to point to instances owned in struct
       locus_sampling. */
    struct dir_points *dist[2];
};


struct alpha_pair {
    unsigned b1, b2;
};

void
print_cache_stats();

/* initializes local 'cache' variable, and calls the init functions
   for services that it depends on. */
void
dirichlet_diff_cache_init(struct dirichlet_diff_params dd_par,
                          struct binomial_est_params be_par,
                          struct dir_cache_params dc_par,
                          struct bam_filter_params bf_par,
                          struct bam_scanner_info *reader_buf,
                          unsigned n_max_reading,
                          unsigned long max_input_mem,
                          unsigned n_threads);

void
dirichlet_diff_cache_free();


void
dir_diff_cache_thread_init();


void
dir_diff_cache_thread_free();


// void prepopulate_bounds_keys(unsigned n_threads);

unsigned
find_cacheable_permutation(const unsigned *a, const unsigned *b, 
                           const unsigned *lim, unsigned *permutation);


/* For two dirichlet distributions A = { x+p, a2+p, p, p } and B = {
   b1+p, b2+p, p, p }, virtually find the values of x in [0, max1)
   that denote the intervals where A and B are UNCHANGED, and where
   they are AMBIGUOUS.  The expected underlying pattern as a function
   of x is some number of C, AC, A, AU, U, AU, A, AC, C.  (see enum
   fuzzy_state in binomial_est.h) */
void
initialize_est_bounds(unsigned a2, unsigned b1, unsigned b2,
                      struct bound_search_params *bpar,
                      struct binomial_est_bounds *beb);


/* test two dirichlets based on their counts. Use a thread-safe
   caching mechanism to update the cache as necessary. sets cache_hit
   if the cache was used. */
enum fuzzy_state
cached_dirichlet_diff(unsigned *alpha1, unsigned *alpha2,
                      struct bound_search_params *bsp,
                      unsigned *cache_hit);


#endif /* _DIR_DIFF_CACHE_H */
