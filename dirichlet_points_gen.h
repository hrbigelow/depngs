#ifndef _DIRICHLET_POINTS_GEN_H
#define _DIRICHLET_POINTS_GEN_H

#include "defs.h"
#include <gsl/gsl_rng.h>

#include "binomial_est.h"

/* number of points or weights generated at a time */
#define GEN_POINTS_BATCH 32

struct dir_points_par
{
    double alpha[NUM_NUCS];
    gsl_rng *randgen;
};

/* Generate GEN_POINTS_BATCH points using par to parameterize the
   distribution */
void gen_dirichlet_points_wrapper(const void *par, POINT *points);

/* generate GEN_POINTS_BATCH 'reference' points, representing the
   corner of the simplex corresponding to the reference base, or a
   point outside the simplex for reference 'N'.  This external point
   will serve as an 'always different' point. */
void gen_reference_points_wrapper(const void *par, POINT *points);


struct calc_post_to_dir_par
{
    struct packed_counts *post_counts;
    double proposal_alpha[NUM_NUCS];
};

/* Generate GEN_POINTS_BATCH weights (ratio of posterior to dirichlet) */
void calc_post_to_dir_ratio(POINT *points, const void *par, double *weights);

/* Generate GEN_POINTS_BATCH dummy weights of value 1 */
void calc_dummy_ratio(POINT * /* unused */, const void * /* unused */,
                      double *weights);


#endif /* _DIRICHLET_POINTS_GEN_H */
