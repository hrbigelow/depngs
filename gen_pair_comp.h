#ifndef _GEN_PAIR_COMP_H
#define _GEN_PAIR_COMP_H

#include <gsl/gsl_rng.h>

struct pair_comp {
    double c1[4];
    double c2[4];
};


/* generate a pair of points in pair_comp, in barycentric coordinates
   that are dist apart in the regular unit simplex coordinates. the
   first point is sampled from a Dirichlet with alpha. */
struct pair_comp gen_pair_comp(double *alpha, double dist, gsl_rng *rg);


#endif /* _GEN_PAIR_COMP_H */
