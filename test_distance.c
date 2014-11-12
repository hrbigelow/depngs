#include <gsl/gsl_rng.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_randist.h>
#include <math.h>
#include <time.h>
#include <string.h>

#include "spatial_search.h"

/*
  In 1, 2, or 3 dimensions, simulate N points in the unit hypercube
  using a uniform distribution.  Then, calculate the average
  nearest-neighbor distance.

 */

#define MIN(a, b) ((a) < (b) ? (a) : (b))

struct wpoint {
    double *x;
    double wgt;
};

int double_comp(const void *a, const void *b)
{
    double av = *(double *)a, bv = *(double *)b;
    return av < bv ? -1 : ((av == bv) ? 0 : 1);
}


static unsigned cmp_dim = 0;

// compare wpoint *pa and *pb
int pcmp(const void *pa, const void *pb)
{
    const struct wpoint 
        *a = (struct wpoint *)pa, 
        *b = (struct wpoint *)pb;

    return 
        a->x[cmp_dim] < b->x[cmp_dim] ? -1 : a->x[cmp_dim] > b->x[cmp_dim] ? 1 : 0;
}


int wcmp(const void *pa, const void *pb)
{
    const struct wpoint 
        *a = (struct wpoint *)pa, 
        *b = (struct wpoint *)pb;

    return a->wgt < b->wgt ? -1 : a->wgt > b->wgt ? 1 : 0;
}

#define SQRT2 1.414213562373095048801688724209
#define SQRT3 1.732050807568877293527446341505

static double simplex2[] = {
    -1.0 / 2.0,
    1.0 / 2.0
};

static double simplex3[] = {
    1.0 / 2.0, 0,
    -1.0 / 2.0, 0,
    0, SQRT3 / 2.0
};

static double simplex4[] = {

    1,  0, -1.0 / SQRT2,
    -1, 0, -1.0 / SQRT2,
    0,  1,  1.0 / SQRT2,
    0, -1,  1.0 / SQRT2
};


static double LO_BOUND = 0.01;
static double HI_BOUND = 0.99;

/* N normalized components fit into N-1 dimensions.  The simplex
   corners are N points, each of N-1 dimensions.  The result is an N-1
   dimensional 'point' */
void comp_to_simplex(double *comp, double *point, int ncomp)
{
    int ci, di, ndim = ncomp - 1;
    double *p;
    switch (ncomp)
    {
    case 2: p = simplex2; break;
    case 3: p = simplex3; break;
    case 4: p = simplex4; break;
    }

    /* iterate over the N simplex corners.  For each one, add its
       contribution to each of the N-1 dimensions */
    for (di = 0; di != ndim; ++di)
    {
        point[di] = 0;
        for (ci = 0; ci != ncomp; ++ci)
            point[di] += p[ci * ndim + di] * comp[ci];
    }
}


/*
double euclidean_dist(double *p1, double *p2, int ncomp, int *ndim)
{
    int d;
    double sq_dist = 0;
    *ndim = ncomp;

    for (d = 0; d != ncomp; ++d)
        sq_dist += gsl_pow_2(p1[d] - p2[d]);

    return sqrt(sq_dist);
}
*/

static double *global_alpha;

double dirichlet_pdf(double *point, int ndim)
{
    return gsl_ran_dirichlet_pdf(ndim, global_alpha, point);
}


/* points are a set of weighted points sampled from some
   non-normalized distribution P(*), each weighted by the value of
   P(*) at that point.  For each point p in points, estimate its
   'neighborhood volume' v by finding its G nearest neighbors.
   Calculate mass as v * P(p) for each point.  Finally, return median
   and mean of masses filtered somehow... */
void partite_mass(struct wpoint *points, unsigned npoints, unsigned G,
                  double min_quantile, double max_quantile,
                  unsigned sim)
{
    /* build the marked_point arrays */
    struct marked_point **mpoints[NDIM];

    struct marked_point 
        *mbuf = (struct marked_point *)malloc(npoints * sizeof(struct marked_point)),
        *mcur,
        *mend = mbuf + npoints;

    /* initialize the content from wpoints to marked_points */
    struct wpoint *wp = points;
    mcur = mbuf;
    while (mcur != mend)
    {
        memcpy(mcur->p, wp->x, sizeof(mcur->p));
        mcur->dist_val = wp->wgt;
        mcur->center = NULL; /* this is the only thing we need to mark this as 'uninitialized' */
        ++mcur;
        ++wp;
    }

    unsigned p, d;

    /* initialize all pointers in mpoints[] */
    for (d = 0; d != NDIM; ++d)
    {
        mpoints[d] = (struct marked_point **)malloc(npoints * sizeof(struct marked_point *));
        struct marked_point
            **mpcur = mpoints[d],
            **mpend = mpcur + npoints;
        mcur = mbuf;
        while (mpcur != mpend)
            *mpcur++ = mcur++;
    }

    /* sort the points by individual dimension */
    for (d = 0; d != NDIM; ++d)
    {
        set_cmp_dim(d);
        qsort(mpoints[d], npoints, sizeof(struct marked_point *), marked_point_comp);
    }    

    /* for each point, do the spatial search, and update the metrics
       of interest */
    struct marked_point *head;
    unsigned total_used = 0;
    double *volume = (double *)malloc(npoints * sizeof(double *));
    double *weight = (double *)malloc(npoints * sizeof(double *));
    double *mass = (double *)malloc(npoints * sizeof(double *));

    for (p = 0; p != npoints; ++p)
    {
        head = spatial_search(mpoints, npoints, mpoints[0][p], G);
        /* estimate the neighborhood volume somehow */
        if (! head)
            continue;

        volume[total_used] = estimate_volume(head, G);
        weight[total_used] = head->center->dist_val;
        mass[total_used] = volume[total_used] * weight[total_used];
        ++total_used;
    }

    /* report raw metrics */
    if (sim == 0)
        fprintf(stdout, "%s\t%s\t%s\t%s\t%s\t%s\n",
                "category",
                "num_points",
                "sim#",
                "volume",
                "weight",
                "mass");
    
    for (p = 0; p != total_used; ++p)
    {
        fprintf(stdout, "%u\t%u\t%u\t%g\t%g\t%g\n",
                p,
                total_used,
                sim,
                volume[p],
                weight[p],
                mass[p]);
    }

    /* calculate metrics of interest */
    qsort(volume, total_used, sizeof(double), double_comp);
    qsort(weight, total_used, sizeof(double), double_comp);
    qsort(mass, total_used, sizeof(double), double_comp);

    double medians[] = { 
        volume[total_used / 2],
        weight[total_used / 2],
        mass[total_used / 2]
    };
    
    unsigned
        min_p = (unsigned)floor(total_used * min_quantile),
        max_p = (unsigned)floor(total_used * max_quantile) - 1,
        num_used_metrics = max_p - min_p;
        
    double means[] = { 0, 0, 0 };
    for (p = min_p; p != max_p; ++p)
    {
        means[0] += volume[p];
        means[1] += weight[p];
        means[2] += mass[p];
    }
    means[0] /= (double)num_used_metrics;
    means[1] /= (double)num_used_metrics;
    means[2] /= (double)num_used_metrics;
    
    fprintf(stdout, "%s\t%u\t%u\t%g\t%g\t%g\n", "mean", 
            total_used, sim, means[0], means[1], means[2]);

    fprintf(stdout, "%s\t%u\t%u\t%g\t%g\t%g\n", "median", 
            total_used, sim, medians[0], medians[1], medians[2]);

    /* clean up */
    for (d = 0; d != NDIM; ++d)
        free(mpoints[d]);

    free(mbuf);
    free(volume);
    free(weight);
    free(mass);
}



/* compute expected value of nearest neighbor distance, using all
   points inside the 8th square (square with 1/2 of the side length of
   its containing square, centered in the middle) */
double mean_nn_dist(struct wpoint *points, int npoints,
                    int ncomp, 
                    double (*dist_fcn)(double *p1, double *p2, int ncomp, int *ndim),
                    double *median_dist, 
                    double *median_wgt,
                    int *npoints_used)
{
    double *dists = (double *)malloc(npoints * sizeof(double));

    /* initialize the points */
    struct wpoint *pl, *pr, *pp, *pbeg = points, *pend = pbeg + npoints;

    int d;
    int ndim;

    /* find the ranges of weight that are safe */
    qsort(points, npoints, sizeof(struct wpoint), wcmp);
    double lo_weight = points[(int)floor(npoints * LO_BOUND)].wgt;
    double hi_weight = points[(int)floor(npoints * HI_BOUND)].wgt;
    
    cmp_dim = 0;
    qsort(points, npoints, sizeof(struct wpoint), pcmp);
    /* struct wpoint **lb = lower_bound(ppoints, npoints, sizeof(double *), &p_lo, pcmp); */
    /* struct wpoint **ub = upper_bound(ppoints, npoints, sizeof(double *), &p_hi, pcmp); */

    int left_more = 1, right_more = 1, total_count = 0, n = 0;
    double nn_dist, cur_dist, total_dist = 0, total_prod = 0;

    for (pp = pbeg; pp != pend; ++pp)
    {
        pl = pp;
        pr = pp;
        int out_of_bounds = 0;
        for (d = 0; d != ncomp; ++d)
            if (pp->wgt < lo_weight || pp->wgt > hi_weight)
            {
                out_of_bounds = 1;
                break;
            }

        /* if (out_of_bounds) continue; */

        left_more = pp != pbeg;
        right_more = pp != pend;
        nn_dist = DBL_MAX;

        while (left_more || right_more)
        {
            if (left_more)
            {
                if (--pl == pbeg || fabs(pl->x[cmp_dim] - pp->x[cmp_dim]) > nn_dist)
                    left_more = 0;

                else
                    nn_dist = ((cur_dist = dist_fcn(pl->x, pp->x, ncomp, &ndim)) < nn_dist)
                        ? cur_dist : nn_dist;
            }
            if (right_more)
            {
                if (++pr == pend || fabs(pr->x[cmp_dim] - pp->x[cmp_dim]) > nn_dist)
                    right_more = 0;
                
                else
                    nn_dist = ((cur_dist = dist_fcn(pr->x, pp->x, ncomp, &ndim)) < nn_dist)
                        ? cur_dist : nn_dist;
            }
        }

        dists[n++] = nn_dist;
        if (! out_of_bounds && pl != pbeg && pr != pend)
        {
            total_dist += nn_dist;
            total_prod += pow(nn_dist, ndim) * pp->wgt;
            total_count++;
        }
        if ((pend - pp) % 10000 == 0)
            fprintf(stdout, "%lu points left.\n", pend - pp);
    }
    qsort(dists, npoints, sizeof(double), double_comp);
    *median_dist = dists[npoints / 2];

    qsort(points, npoints, sizeof(struct wpoint), wcmp);
    *median_wgt = points[npoints / 2].wgt;

    *npoints_used = total_count;

    free(dists);
    return total_prod / (double)total_count;
}


void die(const char *msg, int exitval)
{
    fprintf(stdout, "%s\n", msg);
    exit(exitval);
}

int main(int argc, char **argv)
{
    int npoints = atoi(argv[1]);
    unsigned G = atoi(argv[2]);
    unsigned nsim = atoi(argv[3]);
    char *str_alpha = argv[4];

    unsigned ncomp = NDIM + 1;

    global_alpha = (double *)malloc(npoints * 4 * sizeof(double));

    sscanf(str_alpha, "%lf,%lf,%lf,%lf", &global_alpha[0], &global_alpha[1], &global_alpha[2], &global_alpha[3]);

    int n_dirichlet_dim = NDIM + 1;

    /* double total_volume = 1; */
    /* double points_per_unit_volume = (double)npoints / total_volume; */

    gsl_rng *rand = gsl_rng_alloc(gsl_rng_taus);
    struct timespec now;
    clock_gettime(CLOCK_REALTIME, &now);
    gsl_rng_set(rand, now.tv_sec);

    struct wpoint *points = (struct wpoint *)malloc(npoints * sizeof(struct wpoint));
    double *pb, *buf = (double *)malloc(npoints * NDIM * sizeof(double));

    if (! points)
        die("Couldn't allocate that many points.", 1);

    struct wpoint *pp, *pe = points + npoints;
    /* double mean_vol_times_wgt = 0, median_dist, median_wgt; */
    
    for (pp = points, pb = buf; pp != pe; ++pp, pb += NDIM)
        pp->x = pb;

    /* alternative simulation of points using a Dirichlet */
    double *comp = (double *)malloc(ncomp * sizeof(double));

    /* do X simulations */
    unsigned sim;
    for (sim = 0; sim != nsim; ++sim)
    {
        for (pp = points; pp != pe; ++pp)
        {
            gsl_ran_dirichlet(rand, n_dirichlet_dim, global_alpha, comp);
            pp->wgt = dirichlet_pdf(comp, ncomp);
            comp_to_simplex(comp, pp->x, ncomp);
        }
        
        partite_mass(points, npoints, G, LO_BOUND, HI_BOUND, sim);
    }
    /*
    int npoints_used;
    mean_vol_times_wgt = mean_nn_dist(points, npoints, NDIM, 
                                      euclidean_dist, &median_dist, 
                                      &median_wgt, &npoints_used);

    fprintf(stdout, 
            "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n",
            "ncomp",
            "NDIM",
            "npoints",
            "npoints_used",
            "mean",
            "mean_times_density",
            "median_dist",
            "median_wgt",
            "md_mw_times_density");

    fprintf(stdout, 
            "%i\t%i\t%i\t%i\t%g\t%g\t%g\t%g\t%g\n",
            ncomp, NDIM, npoints, npoints_used, 
            mean_vol_times_wgt,
            mean_vol_times_wgt * points_per_unit_volume, 
            median_dist,
            median_wgt,
            pow(median_dist, NDIM) * median_wgt * points_per_unit_volume);
    */
    

    free(points);
    free(global_alpha);
    free(buf);
    free(comp);
    gsl_rng_free(rand);

    return 0;
}
