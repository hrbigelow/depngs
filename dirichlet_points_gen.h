#ifndef _DIRICHLET_POINTS_GEN_H
#define _DIRICHLET_POINTS_GEN_H

#include "defs.h"
#include <gsl/gsl_rng.h>

#include "batch_pileup.h"

/* number of points or weights generated at a time */
#define GEN_POINTS_BATCH 32

struct points_gen_par
{
    unsigned alpha_counts[NUM_NUCS];
    unsigned alpha_perm[NUM_NUCS]; /* the permutation that was applied to get alpha_counts */
    gsl_rng *randgen;
    struct bqs_count *observed;
    unsigned n_observed;
    unsigned min_base_quality;
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


void alloc_distrib_points(struct distrib_points *dpts);

void free_distrib_points(struct distrib_points *dpts);


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


/* initializes error_probability and alpha_prior.  (no allocation
   needed) */
void
dirichlet_points_gen_init(double _alpha_prior);

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


/* Generate GEN_POINTS_BATCH dummy weights of value 1 */
void
calc_dummy_logratio(POINT * /* unused */, const void * /* unused */,
                    double *weights);

/* exponentiate vals, scaling to avoid underflow or overflow */
void
batch_scaled_exponentiate(double *val, unsigned n_val);



#endif /* _DIRICHLET_POINTS_GEN_H */
