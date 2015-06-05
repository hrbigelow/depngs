#include "gen_pair_comp.h"
#include "simplex.h"

#include <math.h>
#include <string.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_math.h>

/* Functions to generate a pair of composition points (in barycentric
   coordinates).  Generate the first one from a dirichlet with alpha.
   Before attempting to generate a second one, evaluate what the
   maximal distance is, which is the distance of  */


/* These are used just for the purpose of simulating dist == 1
   pairs. */
static double bary_corners[][4] = {
    { 1, 0, 0, 0 },
    { 0, 1, 0, 0 },
    { 0, 0, 1, 0 },
    { 0, 0, 0, 1 }
};


void unit_sphere3(double *x, gsl_rng *rg)
{
    double n[3];
    unsigned i;
    for (i = 0; i != 3; ++i) n[i] = gsl_ran_gaussian(rg, 1.0);
    double denom = sqrt(n[0] * n[0] + n[1] * n[1] + n[2] * n[2]);
    for (i = 0; i != 3; ++i) x[i] = n[i] / denom;
}


double dist4(double *a, double *b)
{
    return sqrt(gsl_pow_2(a[0] - b[0])
                + gsl_pow_2(a[1] - b[1])
                + gsl_pow_2(a[2] - b[2])
                + gsl_pow_2(a[3] - b[3]));
}

#define ONE_OVER_SQRT2 0.70710678118654752440

double dist4_scaled(double *a, double *b)
{
    return dist4(a, b) * ONE_OVER_SQRT2;
}

/* return 1 if the point is a valid barycentric point (normalized
   positive) */
unsigned barycentric(double *p)
{
    double sum = p[0] + p[1] + p[2] + p[3];
    return (sum - 1.0 < 1e-10 && p[0] >= 0 && p[1] >= 0 && p[2] >= 0 && p[3] >= 0);
}

struct pair_comp gen_pair_comp(double *alpha, double dist, gsl_rng *rg)
{
    struct pair_comp pc;
    double s1[3], min_v, max_dist = -1.0;
    unsigned i, min_i;
    double alpha_low[] = { 0.1, 0.1, 0.1, 0.1 };

    if (dist == 1)
    {
        unsigned long c1 = gsl_rng_uniform_int(rg, 4);
        unsigned long c2 = (c1 + 1 + gsl_rng_uniform_int(rg, 3)) % 4;
        memcpy(pc.c1, bary_corners[c1], sizeof(pc.c1));
        memcpy(pc.c2, bary_corners[c2], sizeof(pc.c2));
    }
    else if (dist == 0)
    {
        gsl_ran_dirichlet(rg, 4, (dist > 0.9 ? alpha_low : alpha), pc.c1);
        memcpy(pc.c2, pc.c1, sizeof(pc.c2));
    }
    else if (dist > 0.95)
    {
        fprintf(stderr, "Cannot handle distances not equal to 1 but > 0.95\n");
        exit(1);
    }
    else
    {
        /* generate a valid c1 */
        double s_corner[3];
        while (max_dist < dist * (1.0 / 0.95))
        {
            gsl_ran_dirichlet(rg, 4, (dist > 0.9 ? alpha_low : alpha), pc.c1);
            min_v = pc.c1[0], min_i = 0;
            for (i = 1; i != 4; ++i)
                if (min_v > pc.c1[i])
                    min_v = pc.c1[i], min_i = i;
            barycentric_to_simplex(pc.c1, s1);
            barycentric_to_simplex(bary_corners[min_i], s_corner);
            max_dist = dist3(s_corner, s1);
        }
    
        /* sample from the unit sphere until we find a point in the
           simplex */
        double sph[3], s1p[3];

        while (1)
        {
            unit_sphere3(sph, rg);
            for (i = 0; i != 3; ++i) s1p[i] = s1[i] + dist * sph[i];
            if (inside_simplex(s1p)) break;
        }
        simplex_to_barycentric(s1p, pc.c2);
    }
    /* fprintf(stderr,  */
    /*         "%5.3f,%5.3f,%5.3f,%5.3f\t" */
    /*         "%5.3f,%5.3f,%5.3f,%5.3f\t" */
    /*         "%f\t%f\n",  */
    /*         pc.c1[0], pc.c1[1], pc.c1[2], pc.c1[3], */
    /*         pc.c2[0], pc.c2[1], pc.c2[2], pc.c2[3], */
    /*         dist,  */
    /*         dist4_scaled(pc.c1, pc.c2) */
    /*         ); */
    return pc;
}
