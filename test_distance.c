#include <gsl/gsl_rng.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_randist.h>
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

static unsigned cmp_dim = 0;

/* qsort presents the comparison function with the addresses of two
   elements.*/
int pcmp(const void *pa, const void *pb)
{
    const double *a = *(const double **)pa, *b = *(const double **)pb;
    return 
        a[cmp_dim] < b[cmp_dim] ? -1 : a[cmp_dim] > b[cmp_dim] ? 1 : 0;
}


/* returns the left-most pointer p in [base, base + nmemb * size)
   where compar(p, val) >= 0 */
void *lower_bound(void *base, size_t nmemb, size_t size, 
                  void *val, /* address of the value to be fed into the comparison function */
                  int(*compar)(const void *, const void *))
{
    size_t half;
    void *first = base, *middle;
    while (nmemb > 0)
    {
        half = nmemb >> 1;
        middle = first;
        middle += half * size;
        if (compar(middle, val) < 0)
        {
            first = middle;
            first += size;
            nmemb = nmemb - half - 1;
        }
        else
            nmemb = half;
    }
    return first;
}


/* returns the right-most pointer p in [base, base + nmemb * size)
   where compar(val, p) < 0 or p == base + nmemb * size if this is not
   true of any element. */
void *upper_bound(void *base, size_t nmemb, size_t size, 
                  void *val, /* address of the value to be fed into the comparison function */
                  int(*compar)(const void *, const void *))
{
    size_t half;
    void *first = base, *middle;
    while (nmemb > 0)
    {
        half = nmemb >> 1;
        middle = first;
        middle += half * size;
        if (compar(val, middle) < 0)
            nmemb = half;
        else
        {
            first = middle;
            first += size;
            nmemb = nmemb - half - 1;
        }
    }
    return first;
}


double euclidean_dist(double *p1, double *p2, size_t ndim)
{
    int d;
    double sq_dist = 0;
    for (d = 0; d != ndim; ++d)
        sq_dist += gsl_pow_2(p1[d] - p2[d]);

    return sqrt(sq_dist);
}


/* compute expected value of nearest neighbor distance, using all
   points inside the 8th square (square with 1/2 of the side length of
   its containing square, centered in the middle) */
double mean_nn_dist(double *points, int npoints, double half_side, int ndim)
{
    double **ppoints = (double **)malloc(npoints * sizeof(double *));
    /* initialize the points */
    double **pl, **pr;
    double *p = points;
    double **pp = ppoints;
    while (pp != ppoints + npoints)
    {
        *pp++ = p;
        p += ndim;
    }
    cmp_dim = 0;
    qsort(ppoints, npoints, sizeof(double *), pcmp);
    double lo_val = -half_side / 2.0, hi_val = half_side / 2.0;
    double *p_lo = &lo_val, *p_hi = &hi_val;

    double **lb = lower_bound(ppoints, npoints, sizeof(double *), &p_lo, pcmp);
    double **ub = upper_bound(ppoints, npoints, sizeof(double *), &p_hi, pcmp);

    int left_more = 1, right_more = 1, total_count = 0;
    double nn_dist, cur_dist, total_dist = 0;
    for (pp = lb; pp != ub; ++pp)
    {
        pl = pp;
        pr = pp;
        left_more = pp != ppoints;
        right_more = pp != ppoints + npoints;
        nn_dist = DBL_MAX;

        while (left_more || right_more)
        {
            if (left_more)
            {
                if (--pl == ppoints || fabs((*pl)[cmp_dim] - (*pp)[cmp_dim]) > nn_dist)
                    left_more = 0;

                else
                    nn_dist = ((cur_dist = euclidean_dist(*pl, *pp, ndim)) < nn_dist)
                        ? cur_dist : nn_dist;
            }
            if (right_more)
            {
                if (++pr == ppoints + npoints || fabs((*pr)[cmp_dim] - (*pp)[cmp_dim]) > nn_dist)
                    right_more = 0;
                
                else
                    nn_dist = ((cur_dist = euclidean_dist(*pr, *pp, ndim)) < nn_dist)
                        ? cur_dist : nn_dist;
            }
        }
        if (pl != ppoints && pr != ppoints + npoints)
        {
            total_dist += nn_dist;
            ++total_count;
        }
        if ((ub - pp) % 10000 == 0)
            fprintf(stdout, "%i points left.\n", ub - pp);
    }
    free(ppoints);
    return total_dist / (double)total_count;
}


void die(char *msg, int exitval)
{
    fprintf(stdout, "%s\n", msg);
    exit(exitval);
}

int main(int argc, char **argv)
{
    /* */
    int ndim = atoi(argv[1]);
    int points_per_unit_vol = atoi(argv[2]);
    double half_side = atof(argv[3]);
    int npoints = points_per_unit_vol * gsl_pow_int((int)floor(half_side) * 2.0, ndim);

    gsl_rng *rand = gsl_rng_alloc(gsl_rng_taus);

    int c, p, d;
    double *points = (double *)malloc(npoints * ndim * sizeof(double));
    if (! points)
        die("Couldn't allocate that many points.", 1);

    double *pp, *pe = points + npoints * ndim;
    /* double *comps = (double *)malloc(ncomps * sizeof(double)); */
    double mean = 0;

    for (pp = points; pp != pe; pp += ndim)
        for (d = 0; d != ndim; ++d)
            pp[d] = gsl_ran_flat(rand, -half_side, half_side);

    mean = mean_nn_dist(points, npoints, half_side, ndim);

    double density = pow(points_per_unit_vol, 1.0 / (double)ndim);
    fprintf(stdout, "ndim\tnpoints\tncomps\tmedian\tmean\tmean_times_density\n");
    fprintf(stdout, "%i\t%i\t%g\t%g\n", ndim, npoints, mean, mean * density);

    free(points);
    /* free(comps); */
    gsl_rng_free(rand);

    return 0;
}
