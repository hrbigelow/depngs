#ifndef _DIRICHLET_DIFF_CACHE_H
#define _DIRICHLET_DIFF_CACHE_H

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
    struct distrib_points *dist[2];
};


extern unsigned alpha_packed_limits[];


/* struct alpha_packed_large { */
/*     unsigned a0 :24; /\* 16,777,216 *\/ */
/*     unsigned a1 :20; /\*  1,048,576 *\/ */
/*     unsigned a2 :12; /\*      4,096 *\/ */
/*     unsigned a3 :8;  /\*        256 *\/ */
/* }; */


khint64_t
pack_alpha64(unsigned a0, unsigned a1, unsigned a2, unsigned a3);


void
unpack_alpha64(khint64_t k, unsigned *c);


khint64_t
pack_bounds(unsigned a2, unsigned b1, unsigned b2);


void
unpack_bounds(khint64_t k, unsigned *b);

/* layout is:  
   c[0]: a0:24, a3:8. 
   c[1]: a1:20, a2:12 
*/
/* union alpha_large_key { */
/*     uint32_t c[2]; */
/*     uint64_t key; */
/* }; */


/* union bounds_key { */
/*     struct { */
/*         unsigned a2:20; /\*  1,048,576 *\/ */
/*         unsigned b1:24; /\* 16,777,216 *\/ */
/*         unsigned b2:20; /\*  1,048,576 *\/ */
/*     } f; */
/*     int64_t key; */
/* }; */



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
                          void **reader_pars,
                          struct contig_region *qbeg,
                          struct contig_region *qend,
                          unsigned n_max_reading,
                          unsigned long max_input_mem,
                          unsigned n_threads);

void
dirichlet_diff_cache_free();

// void prepopulate_bounds_keys(unsigned n_threads);

void
find_cacheable_permutation(const unsigned *a, const unsigned *b, 
                           const unsigned *lim,
                           unsigned *permutation, 
                           unsigned *perm_found);

void
update_points_gen_params(struct distrib_points *dpts,
                         unsigned *alpha_counts,
                         unsigned *permutation);


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
cached_dirichlet_diff(unsigned *a_counts,
                      unsigned *b_counts,
                      struct bound_search_params *bpar,
                      unsigned *cacheable,
                      unsigned *cache_was_set);


#endif /* _DIRICHLET_DIFF_CACHE_H */
