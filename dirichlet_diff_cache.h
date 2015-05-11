#ifndef _DIRICHLET_DIFF_CACHE_H
#define _DIRICHLET_DIFF_CACHE_H

#include "binomial_est.h"

#include <stdio.h>

/* this is instantiated inside locus_sampling (one instance per
   (thread x sample))*/
struct distrib_points {
    struct points_gen pgen;
    struct points_buf points;
    struct weights_buf weights;
};


void alloc_distrib_points(struct distrib_points *dpts);

void free_distrib_points(struct distrib_points *dpts);


/* one instance per thread. */
struct binomial_est_params {
    enum fuzzy_state query_state;
    unsigned use_low_beta;
    double query_beta;
    unsigned points_hash_frozen;
    struct distrib_points *dist[2];
};


struct alpha_pair {
    unsigned b1, b2;
};

void print_primary_cache_size();

unsigned points_hash_frozen();

void print_cache_stats();

void dirichlet_diff_init(unsigned max1, unsigned max2, 
                         unsigned batch_size,
                         double post_confidence,
                         double beta_confidence,
                         double min_dirichlet_dist,
                         unsigned max_sample_points,
                         unsigned long max_dir_cache_size,
                         unsigned long max_secondary_cache_size,
                         unsigned n_threads);
void dirichlet_diff_free();

void find_cacheable_permutation(const unsigned *a, const unsigned *b, 
                                const unsigned *lim,
                                unsigned *permutation, 
                                unsigned *perm_found);

void update_points_gen_params(struct distrib_points *dpts,
                              unsigned *alpha_counts,
                              unsigned *permutation);

/* test two dirichlets based on their counts. Use a thread-safe
   caching mechanism to update the cache as necessary. sets cache_hit
   if the cache was used. */
enum fuzzy_state
cached_dirichlet_diff(unsigned *a_counts,
                      unsigned *b_counts,
                      struct binomial_est_params *bpar,
                      unsigned points_hash_frozen,
                      unsigned *cacheable,
                      unsigned *cache_was_set);


#endif /* _DIRICHLET_DIFF_CACHE_H */


#if 0
void print_beb_bounds(struct binomial_est_params *bpar);

void print_bounds(struct binomial_est_params *bpar);


/* parse diststats header, initializing fields of pset and max1 and max2  */
void parse_diststats_header(FILE *diststats_fh, double *prior_alpha);

void write_diststats_header(FILE *diststats_fh);

void write_diststats_line(FILE *fh,
                          unsigned a2,
                          unsigned b1,
                          unsigned b2,
                          struct binomial_est_bounds *beb);


/* initialize internal bounds_cache */
void parse_diststats_body(FILE *diststats_fh, unsigned max1, unsigned max2);
#endif
