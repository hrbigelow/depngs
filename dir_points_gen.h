#ifndef _DIR_POINTS_GEN_H
#define _DIR_POINTS_GEN_H

#include "defs.h"
#include <gsl/gsl_rng.h>

#include "batch_pileup.h"

/* number of points or weights generated at a time */
#define GEN_POINTS_BATCH 32

struct dirichlet_points_gen_params {
    unsigned min_base_quality;
    unsigned max_sample_points;
    double alpha_prior;
};

/* initializes error_probability and alpha_prior.  (no allocation
   needed) */
void
dirichlet_points_gen_init(struct dirichlet_points_gen_params pg_par);


void
dirichlet_points_gen_free();


void
dir_points_thread_init();


void
dir_points_thread_free();



/* points and weights are allocated to max_sample_points.  n_points
   and n_weights indicate the number of points or weights that are
   current w.r.t alpha. when perm_alpha is updated, n_points and n_weights
   are set to zero to signal that the existing data are out of
   date. */
struct dir_points {
    unsigned perm_alpha[4], perm[4];
    POINT *points_buf, *data;
    size_t n_points;

    struct bqs_count *bqs_ct;
    unsigned n_bqs_ct;
    double *weights;
    size_t n_weights;
};


/* update dp->alpha and dp->perm.  if a change is detected, updates
   data, n_points, and n_weights as appropriate.  after the call, dp
   may be used for on-demand points generation. */
void
dir_points_update_alpha(const unsigned *alpha,
                        const unsigned *perm,
                        struct dir_points *dp);


/* fully populate the points buffer if necessary */
void
dir_points_fill(struct dir_points *dp);


/* */
void
dir_weights_update_terms(struct bqs_count *bqs_ct, unsigned n_bqs_ct,
                         struct dir_points *dp);

/* fully populate the weights buffer if necessary */
void
dir_weights_fill(struct dir_points *dp);


/* caches locus-specific summary data for an individual sample, so
   that it can be re-used in multiple pairings */
struct locus_data {
    unsigned char confirmed_changed;
    struct {
        unsigned char base_ct: 1;
        unsigned char bqs_ct: 1;
        unsigned char indel_ct: 1;
        unsigned char sample_data: 1;
    } init; /* if these flags are set, means the following fields are
               initialized */

    struct dir_points dist;

    struct base_count base_ct;
    struct bqs_count *bqs_ct;
    unsigned n_bqs_ct;
    struct indel_count *indel_ct;
    unsigned n_indel_ct;
    struct pileup_data sample_data;
};

void
alloc_locus_data(struct locus_data *ld);

void
free_locus_data(struct locus_data *ld);

void
reset_locus_data(struct locus_data *ld);


double
get_alpha_prior();


/* generate a complete set of dirichlet points according to cts. */
void
gen_dir_points(unsigned *cts, POINT *points, unsigned n_points);


/* re-order point components to the default permutation { 0, 1, 2, 3
   }.  if dp is using the cache, copy those points to its local
   storage first. */
void
de_permute_points(struct dir_points *dp);


/* Generate GEN_POINTS_BATCH points using par to parameterize the
   distribution */
void
gen_dirichlet_points_wrapper(const void *par, POINT *points);

/* generate GEN_POINTS_BATCH 'reference' points, representing the
   corner of the simplex corresponding to the reference base, or a
   point outside the simplex for reference 'N'.  This external point
   will serve as an 'always different' point. */
void
gen_reference_points_wrapper(const void *par, POINT *points);


/* Generate GEN_POINTS_BATCH weights (ratio of posterior to dirichlet) */
void
calc_post_to_dir_logratio(struct dir_points *dp);


/* populates square_dist_buf with squares of euclidean distances
   between points1 and points2 in the barycentric space (R4,
   normalized positive components).  populates weights_buf with
   product of weights1 and weights2. */
void
compute_wsq_dist(const double *points1, const double *weights1,
                 const double *points2, const double *weights2,
                 size_t n_points,
                 double *square_dist_buf, double *weights_buf);


/* Generate GEN_POINTS_BATCH dummy weights of value 1 */
void
calc_dummy_logratio(POINT * /* unused */, const void * /* unused */,
                    double *weights);

/* exponentiate vals, scaling to avoid underflow or overflow */
void
batch_scaled_exponentiate(double *val, unsigned n_val);



#endif /* _DIR_POINTS_GEN_H */
