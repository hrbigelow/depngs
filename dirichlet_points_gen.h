#ifndef _DIRICHLET_POINTS_GEN_H
#define _DIRICHLET_POINTS_GEN_H

#include "defs.h"
#include <gsl/gsl_rng.h>

#include "binomial_est.h"

/* number of points or weights generated at a time */
#define GEN_POINTS_BATCH 32

struct points_gen_par
{
    unsigned alpha_counts[NUM_NUCS];
    unsigned alpha_perm[NUM_NUCS]; /* the permutation that was applied to get alpha_counts */
    gsl_rng *randgen;
    struct bqs_count *observed;
    unsigned n_observed;
    /* struct packed_counts *post_counts; */
};

void init_dirichlet_points_gen(double _alpha_prior);

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
void calc_post_to_dir_ratio(POINT *points, const void *par, double *weights);

/* Generate GEN_POINTS_BATCH dummy weights of value 1 */
void calc_dummy_ratio(POINT * /* unused */, const void * /* unused */,
                      double *weights);


#endif /* _DIRICHLET_POINTS_GEN_H */
