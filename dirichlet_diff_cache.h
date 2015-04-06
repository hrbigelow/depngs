#ifndef _DIRICHLET_DIFF_CACHE_H
#define _DIRICHLET_DIFF_CACHE_H

#include "binomial_est.h"
#include "metropolis_sampling.h"

struct distrib_points {
    struct points_gen pgen;
    struct points_buf points;
    struct weights_buf weights;
};


void alloc_distrib_points(struct distrib_points *dpts,
                          unsigned max_sample_points);

void free_distrib_points(struct distrib_points *dpts);


/* one instance per thread. */
struct binomial_est_params {
    enum fuzzy_state query_state;
    unsigned use_low_beta;
    double query_beta;
    struct posterior_settings *pset;
    struct distrib_points *dist[2];
    size_t batch_size;
};

enum init_phase { UNSET, PENDING, SET };

/* Describes estimated distance between two Dirichlets with alphas equal to:
   { A1 + p, A2 + p, p, p }
   { B1 + p, B2 + p, p, p }

   Where A1, A2, B1, B2 are integral values.  [unchanged[0],
   unchanged[1]) represents the range of values of A1 where it is
   deemed UNCHANGED.  [ambiguous[0], ambiguous[1]) are the values of
   A1 for which it is deemed AMBIGUOUS (or UNCHANGED, where this
   interval overlaps the unchanged interval).  By construction: 
   0 <= A[0] <= U[0] <= U[1] <= A[1] <= MAX_COUNT1
   
*/
struct binomial_est_bounds {
    enum init_phase state;
    int16_t ambiguous[2];
    int16_t unchanged[2];
};


void dirichlet_diff_init();
void dirichlet_diff_free();


void print_beb_bounds(struct binomial_est_params *bpar);

void print_bounds(struct binomial_est_params *bpar);


/* For two dirichlet distributions A = { x+p, a2+p, p, p } and B = {
   b1+p, b2+p, p, p }, virtually find the values of x in [0,
   MAX_COUNT1) that denote the intervals where A and B are UNCHANGED,
   and where they are AMBIGUOUS.  The expected underlying pattern as a
   function of x is some number of C, AC, A, AU, U, AU, A, AC, C.
   (see enum fuzzy_state in binomial_est.h) */
void initialize_est_bounds(unsigned a2, unsigned b1, unsigned b2,
                           struct binomial_est_params *bpar,
                           struct binomial_est_bounds *beb);


/* test two dirichlets based on their counts. Use a thread-safe
   caching mechanism to update the cache as necessary. */
enum fuzzy_state
cached_dirichlet_diff(unsigned *a_counts,
                      unsigned *b_counts,
                      struct binomial_est_params *bpar);


#endif /* _DIRICHLET_DIFF_CACHE_H */
