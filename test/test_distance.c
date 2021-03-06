#include <gsl/gsl_rng.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_statistics_double.h>
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


double min_component(double *comp, unsigned ncomp)
{
    double c = DBL_MAX, tmp;
    while (ncomp--)
        c = (tmp = *comp++) < c ? tmp : c;
    return c;
}

#define SQRT_THREE_HALVES 1.224744871391589049098642037352

/* the distance of a point in the domain of an N-dimensional
   Dirichlet to the boundary of that domain is determined only by the  */
double dist_to_dirichlet_boundary(double *point, unsigned N)
{
    double dist, mc = min_component(point, N);
    
    switch(N)
    {
    case 2: dist = mc; break;
    case 3: dist = SQRT_THREE_HALVES * mc; break;
    case 4: dist = (2.0 * M_SQRT3 / 3.0) * mc; break;
    }
    return dist;
}

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

static double *global_alpha;

double dirichlet_pdf(double *point, int ndim)
{
    return gsl_ran_dirichlet_pdf(ndim, global_alpha, point);
}


struct metrics
{
    double volume, weight, mass;
};


/* points are a set of weighted points sampled from some
   non-normalized distribution P(*), each weighted by the value of
   P(*) at that point.  For each point p in points, estimate its
   'neighborhood volume' v by finding its G nearest neighbors.
   Calculate mass as v * P(p) for each point.  Finally, return median
   and mean of masses filtered somehow... */
void partite_mass(struct wpoint *points, unsigned npoints, unsigned G,
                  double min_quantile, double max_quantile,
                  unsigned sim,
                  struct metrics *means, /* volume, weight, mass */
                  struct metrics *medians /* volume, weight, mass */)
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

    double dist_to_boundary;
    for (p = 0; p != npoints; ++p)
    {
        dist_to_boundary = dist_to_dirichlet_boundary(mpoints[0][p]->p, NDIM);
        head = spatial_search(mpoints, npoints, mpoints[0][p], dist_to_boundary, G);
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

    medians->volume = volume[total_used / 2];
    medians->weight = weight[total_used / 2];
    medians->mass = mass[total_used / 2];
    
    unsigned
        min_p = (unsigned)floor(total_used * min_quantile),
        max_p = (unsigned)floor(total_used * max_quantile) - 1,
        num_used_metrics = max_p - min_p;
        
    means->volume = 0;
    means->weight = 0;
    means->mass = 0;
    for (p = min_p; p != max_p; ++p)
    {
        means->volume += volume[p];
        means->weight += weight[p];
        means->mass += mass[p];
    }
    means->volume /= (double)num_used_metrics;
    means->weight /= (double)num_used_metrics;
    means->mass /= (double)num_used_metrics;

    fprintf(stdout, "%s\t%u\t%u\t%g\t%g\t%g\n", "mean", 
            total_used, sim, means->volume, means->weight, means->mass);

    fprintf(stdout, "%s\t%u\t%u\t%g\t%g\t%g\n", "median", 
            total_used, sim, medians->volume, medians->weight, medians->mass);

    /* clean up */
    for (d = 0; d != NDIM; ++d)
        free(mpoints[d]);

    free(mbuf);
    free(volume);
    free(weight);
    free(mass);
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

    unsigned ncomp = NDIM;

    global_alpha = (double *)malloc(npoints * 4 * sizeof(double));

    sscanf(str_alpha, "%lf,%lf,%lf,%lf", &global_alpha[0], &global_alpha[1], &global_alpha[2], &global_alpha[3]);

    int n_dirichlet_dim = NDIM;

    gsl_rng *rand = gsl_rng_alloc(gsl_rng_taus);
    struct timespec now;
    clock_gettime(CLOCK_REALTIME, &now);
    gsl_rng_set(rand, now.tv_sec);

    struct wpoint *points = (struct wpoint *)malloc(npoints * sizeof(struct wpoint));
    double *pb, *buf = (double *)malloc(npoints * NDIM * sizeof(double));

    if (! points)
        die("Couldn't allocate that many points.", 1);

    struct wpoint *pp, *pe = points + npoints;
    
    for (pp = points, pb = buf; pp != pe; ++pp, pb += NDIM)
        pp->x = pb;

    /* alternative simulation of points using a Dirichlet */
    double *comp = (double *)malloc(ncomp * sizeof(double));

    struct metrics 
        *means = (struct metrics *)malloc(nsim * sizeof(struct metrics)), 
        *pmeans = means,
        *medians = (struct metrics *)malloc(nsim * sizeof(struct metrics)),
        *pmedians = medians;
    
    /* do X simulations */
    unsigned sim;
    for (sim = 0; sim != nsim; ++sim)
    {
        for (pp = points; pp != pe; ++pp)
        {
            gsl_ran_dirichlet(rand, n_dirichlet_dim, global_alpha, pp->x);
            pp->wgt = dirichlet_pdf(pp->x, n_dirichlet_dim);
            /* comp_to_simplex(comp, pp->x, ncomp); */
        }
        
        partite_mass(points, npoints, G, LO_BOUND, HI_BOUND, sim, pmeans++, pmedians++);
    }

    struct metrics sdmeans, sdmedians;
    sdmeans.volume = gsl_stats_sd((double *)means, 3, nsim);
    sdmeans.weight = gsl_stats_sd(((double *)means) + 1, 3, nsim);
    sdmeans.mass = gsl_stats_sd(((double *)means) + 2, 3, nsim);

    sdmedians.volume = gsl_stats_sd((double *)medians, 3, nsim);
    sdmedians.weight = gsl_stats_sd(((double *)medians) + 1, 3, nsim);
    sdmedians.mass = gsl_stats_sd(((double *)medians) + 2, 3, nsim);

    fprintf(stdout, "%s\t%u\t%i\t%g\t%g\t%g\n",
            "means_sd",
            nsim,
            -1,
            sdmeans.volume,
            sdmeans.weight,
            sdmeans.mass);

    fprintf(stdout, "%s\t%u\t%i\t%g\t%g\t%g\n",
            "medians_sd",
            nsim,
            -1,
            sdmedians.volume,
            sdmedians.weight,
            sdmedians.mass);
    
    free(points);
    free(global_alpha);
    free(buf);
    free(comp);
    free(means);
    free(medians);
    gsl_rng_free(rand);

    return 0;
}
