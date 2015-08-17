#ifndef _DIRICHLET_DIFF_CACHE_H
#define _DIRICHLET_DIFF_CACHE_H

#include "binomial_est.h"

#include <stdio.h>

/* global configuration parameters for dirichlet_diff_cache */
struct dirichlet_diff_params {
    unsigned pseudo_depth;
    unsigned batch_size;
    double post_confidence;
    double beta_confidence;
    double min_dirichlet_dist;
    unsigned max_sample_points;
    double prior_alpha;
} g_dd_par;

struct binomial_est_bounds {
    /* enum init_phase state; */
    int32_t ambiguous[2];
    int32_t unchanged[2];
};


/* one instance per thread. */
struct binomial_est_params {
    enum fuzzy_state query_state;
    unsigned use_low_beta;
    double query_beta;
    /* unsigned points_hash_frozen; */
    /* unsigned bounds_hash_frozen; */

    /* these two are set to point to instances owned in struct
       locus_sampling. */
    struct distrib_points *dist[2];
};


static unsigned alpha_packed_limits[] = {
    1<<24, 1<<20, 1<<12, 1<<8
};

struct alpha_packed_large {
    unsigned a0 :24; /* 16,777,216 */
    unsigned a1 :20; /*  1,048,576 */
    unsigned a2 :12; /*      4,096 */
    unsigned a3 :8;  /*        256 */
};

union alpha_large_key {
    struct alpha_packed_large c;
    uint64_t key;
};


union bounds_key {
    struct {
        unsigned a2:20; /*     1,048,576 */
        unsigned b1:32; /* 4,294,967,296 */
        unsigned b2:20; /*     1,048,576 */
    } f;
    int64_t key;
};



struct alpha_pair {
    unsigned b1, b2;
};

#if 0
void print_primary_cache_size();

void set_points_hash_flag(unsigned disable);

unsigned freeze_points_hash();
unsigned freeze_bounds_hash();

/* call this if the thread needs to stop writing to the shared data */    
void inactivate_shared_data(unsigned inactivate_points, 
                            unsigned inactivate_bounds);
#endif

void print_cache_stats();

/* initializes local 'cache' variable, and calls the init functions
   for services that it depends on. */
void
dirichlet_diff_cache_init(unsigned pseudo_depth,
                          unsigned batch_size,
                          double post_confidence,
                          double beta_confidence,
                          double prior_alpha,
                          double min_dirichlet_dist,
                          unsigned max_sample_points,
                          void **reader_pars,
                          struct contig_region *qbeg,
                          struct contig_region *qend,
                          unsigned n_max_reading,
                          unsigned long max_input_mem,
                          unsigned n_bounds,
                          unsigned min_ct_keep_bound,
                          unsigned n_point_sets,
                          unsigned n_threads);

void
dirichlet_diff_cache_free();

// void prepopulate_bounds_keys(unsigned n_threads);

void find_cacheable_permutation(const unsigned *a, const unsigned *b, 
                                const unsigned *lim,
                                unsigned *permutation, 
                                unsigned *perm_found);

void update_points_gen_params(struct distrib_points *dpts,
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
                      struct binomial_est_params *bpar,
                      struct binomial_est_bounds *beb);


/* test two dirichlets based on their counts. Use a thread-safe
   caching mechanism to update the cache as necessary. sets cache_hit
   if the cache was used. */
enum fuzzy_state
cached_dirichlet_diff(unsigned *a_counts,
                      unsigned *b_counts,
                      struct binomial_est_params *bpar,
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
