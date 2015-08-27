#ifndef _DIRICHLET_POINTS_GEN_H
#define _DIRICHLET_POINTS_GEN_H

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

struct points_gen_par
{
    unsigned alpha_counts[NUM_NUCS];
    unsigned alpha_perm[NUM_NUCS]; /* the permutation that was applied to get alpha_counts */
    gsl_rng *randgen;
    struct bqs_count *observed;
    unsigned n_observed;
};

struct points_buf {
    POINT *buf, *p; /* p will point to either buf or a hash entry */
    size_t size;
};

struct weights_buf {
    double *buf;
    size_t size, alloc;
};


struct points_gen
{
    void *points_gen_par;
    void (*gen_point)(const void *par, POINT *points);
    void (*weight)(POINT *points, const void *par,
                   double *weights);
};


/* (one instance per (thread x sample))*/
struct distrib_points {
    struct points_gen pgen;
    struct points_buf points;
    struct weights_buf weights;
};


/* initializes error_probability and alpha_prior.  (no allocation
   needed) */
void
dirichlet_points_gen_init(struct dirichlet_points_gen_params pg_par);


void
dirichlet_points_gen_free();


void alloc_distrib_points(struct distrib_points *dpts);

void free_distrib_points(struct distrib_points *dpts);

/* populate points buffer, either by cache retrieval or
   computation. */
void
dirichlet_refresh_points(struct distrib_points *dpts);


/* populate weights buffer by computation.  requires a populated
   points buffer. */
void
dirichlet_refresh_weights(struct distrib_points *dpts);


/* caches locus-specific summary data for an individual sample, so
   that it can be re-used in multiple pairings */
struct locus_data {
    unsigned char confirmed_changed;
    struct {
        unsigned char distp: 1;
        unsigned char base_ct: 1;
        unsigned char bqs_ct: 1;
        unsigned char indel_ct: 1;
        unsigned char sample_data: 1;
    } init; /* if these flags are set, means the following fields are
               initialized */

    struct distrib_points distp;
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


double get_alpha_prior();

/* Generate GEN_POINTS_BATCH points using par to parameterize the
   distribution */
void gen_dirichlet_points_wrapper(const void *par, POINT *points);

/* generate GEN_POINTS_BATCH 'reference' points, representing the
   corner of the simplex corresponding to the reference base, or a
   point outside the simplex for reference 'N'.  This external point
   will serve as an 'always different' point. */
void gen_reference_points_wrapper(const void *par, POINT *points);


/* Generate GEN_POINTS_BATCH weights (ratio of posterior to dirichlet) */
void
calc_post_to_dir_logratio(POINT *points, const void *par, double *weights);


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



#endif /* _DIRICHLET_POINTS_GEN_H */
