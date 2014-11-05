#include <gsl/gsl_rng.h>
#include <gsl/gsl_math.h>
#include <math.h>


/*
  In 1, 2, or 3 dimensions, simulate N points in the unit hypercube
  using a uniform distribution.  Then, calculate the average
  nearest-neighbor distance.

 */

#define MIN(a, b) ((a) < (b) ? (a) : (b))

int dcomp(const void *a, const void *b)
{
    double av = *(double *)a, bv = *(double *)b;
    return av < bv ? -1 : ((av == bv) ? 0 : 1);
}

int main(int argc, char **argv)
{
    /* */
    int ndim = atoi(argv[1]);
    int npoints = atoi(argv[2]);
    int ncomps = atoi(argv[3]);

    gsl_rng *rand = gsl_rng_alloc (gsl_rng_taus);

    int c, p, d;
    double *points = (double *)malloc(npoints * ndim * sizeof(double));
    double *pp, *pe = points + npoints * ndim;
    double *comps = (double *)malloc(ncomps * sizeof(double));
    double mean = 0;
    c = 0;
    for (c = 0; c != ncomps; ++c)
    {
        for (pp = points; pp != pe; pp += ndim)
            for (d = 0; d != ndim; ++d)
                pp[d] = gsl_rng_uniform(rand);

        for (pk = 0; pk != 10; ++pk)
        {
            /* choose one point at random */
            double *ref_pp = points + gsl_rng_uniform_int(rand, npoints) * ndim;
            
            /* find nearest neighbor distance */
            double nn_dist_squared = 1, dist_squared, diff;
            for (pp = points; pp != pe; pp += ndim)
            {
                dist_squared = 0;
                for (d = 0; d != ndim; ++d)
                {
                    diff = fabs(ref_pp[d] - pp[d]);
                    dist_squared += gsl_pow_2(MIN(diff, 1 - diff));
                }
                nn_dist_squared = ref_pp == pp ? nn_dist_squared : MIN(nn_dist_squared, dist_squared);
            }
            comps[c] = sqrt(nn_dist_squared);
            mean += comps[c];
            ++c;
        }

    mean /= (double)ncomps;

    qsort(comps, ncomps, sizeof(double), dcomp);

    double density = pow(npoints, 1.0 / (double)ndim);
    fprintf(stdout, "ndim\tnpoints\tncomps\tmedian\tmean\tmean_times_density\n");
    fprintf(stdout, "%i\t%i\t%i\t%g\t%g\t%g\n", ndim, npoints, ncomps, comps[ncomps / 2], mean, mean * density);

    free(points);
    free(comps);
    gsl_rng_free(rand);

    return 0;
}
