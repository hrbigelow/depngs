#include <gsl/gsl_randist.h>
#include <gsl/gsl_cdf.h>
#include <stdio.h>
#include <float.h>
#include <string.h>

struct wpoint {
    double x[4], d, o, z; // density and 'other'.
    char dist[8]; // distribution this point came from
    /* z is the normalization constant of this distribution's sample points */
};

struct ival {
    double lo, hi;
};

static unsigned cmp_dim = 0;

int pcmp5(const void *pa, const void *pb)
{
    const struct wpoint *a = *(const struct wpoint **)pa, *b = *(const struct wpoint **)pb;
    return 
        a->x[cmp_dim] < b->x[cmp_dim] 
        ? -1 
        : (a->x[cmp_dim] > b->x[cmp_dim] ? 1 : 0);
}


inline double max(double x, double y) { return x > y ? x : y; }
inline double min(double x, double y) { return x < y ? x : y; }

#define MIN4(a) min(min((a)[0], (a)[1]), min((a)[2], (a)[3]))
#define MAX4(a) max(max((a)[0], (a)[1]), max((a)[2], (a)[3]))

unsigned min4_index(double *a)
{
    unsigned i = 0;
    double m = DBL_MAX, *p = a + 4;
    while (a != p)
    {
        if (*a < m)
        {
            m = *a;
            i = 4 - (p - a);
        }
        ++a;
    }
    return i;
}

unsigned max4_index(double *a)
{
    unsigned i = 0;
    double m = -DBL_MAX, *p = a + 4;
    while (a != p)
    {
        if (*a > m)
        {
            m = *a;
            i = 4 - (p - a);
        }
        ++a;
    }
    return i;
}

// print out the CDF given a number of points
void print_cdf(struct wpoint **points, unsigned npoints,
               struct wpoint *first_point, // needed to give each point an index
               const char *dist_name,
               FILE *cdf_fh)
{
    unsigned dim, i;
    double sum;
    char fmt[] = "%s\t%c\t%20.18le\t%u\t%20.18le\t%20.18le\t%g\t%s\n"; 
    char base[] = "ACGT";
    struct wpoint *s;

    double z = 0;
    for (i = 0; i != npoints; ++i)
        z += points[i]->o / points[i]->z;

    for (dim = 0; dim != 4; ++dim)
    {
        cmp_dim = dim;
        qsort(points, npoints, sizeof(struct wpoints *), pcmp5);
        sum = 0;
        for (i = 0; i != npoints; ++i)
        {
            s = points[i];
                fprintf(cdf_fh, 
                        fmt, dist_name, base[dim], (unsigned)(s - first_point),
                        s->x[dim], (sum += s->o / s->z / z), s->d, s->z, s->dist);
        }
    }

}

#define BETWEEN(x, i) ((x) > (i).lo && (x) < (i).hi)

void ran_dirichlet_clamped(gsl_rng *rand,
                           size_t K, 
                           const struct ival *ival,
                           const double *alpha,
                           double *theta)
{
    size_t i;
    double clamped, v = 0;
    for (i = 0; i != K; ++i)
    {
        clamped = gsl_ran_flat(rand, ival[i].lo, ival[i].hi);
        theta[i] = gsl_cdf_gamma_Pinv(clamped, alpha[i], 1);
        v += theta[i];
    }
    for (i = 0; i != K; ++i)
        theta[i] /= v;
}



int main(int argc, char **argv)
{
    char *str_alpha_p = argv[1], *str_alpha_q = argv[2];
    char *cdf = argv[3], *points = argv[4];

    unsigned npoints_prior = atoi(argv[5]);
    unsigned npoints_post = atoi(argv[6]);
    unsigned npoints_gold = atoi(argv[7]);


    FILE *cdf_fh = fopen(cdf, "w");
    FILE *points_fh = fopen(points, "w");

    double alpha_p[4], alpha_q[4];
    sscanf(str_alpha_p, "%lf,%lf,%lf,%lf", &alpha_p[0], &alpha_p[1], &alpha_p[2], &alpha_p[3]);
    sscanf(str_alpha_q, "%lf,%lf,%lf,%lf", &alpha_q[0], &alpha_q[1], &alpha_q[2], &alpha_q[3]);

    struct wpoint *wsamples = (struct wpoint *)malloc((npoints_prior + npoints_post) * sizeof(struct wpoint));
    struct wpoint **pws = (struct wpoint **)malloc((npoints_prior + npoints_post) * sizeof(struct wpoint *));

    struct wpoint *s, *se1 = wsamples + npoints_prior, *se2 = wsamples + npoints_prior + npoints_post;

    gsl_rng *rand = gsl_rng_alloc(gsl_rng_taus);

    double zp = 0, zq = 0, zr, zri; // zp/zq, zq / zp
    double sum_p = 0, sum_q = 0;

    // populate P(*)
    for (s = wsamples; s != se1; ++s)
    {
        gsl_ran_dirichlet(rand, 4, alpha_p, s->x);
        s->d = gsl_ran_dirichlet_pdf(4, alpha_p, s->x);
        s->o = gsl_ran_dirichlet_pdf(4, alpha_q, s->x);
        strcpy(s->dist, "prior");
        zp += 1.0 / s->d;
    }
    zp /= (double)(se1 - wsamples);

    for (s = wsamples; s != se1; ++s)
        s->z = zp;

    // populate Q(*)
    for (s = se1; s != se2; ++s)
    {
        gsl_ran_dirichlet(rand, 4, alpha_q, s->x);
        s->d = gsl_ran_dirichlet_pdf(4, alpha_q, s->x);
        s->o = gsl_ran_dirichlet_pdf(4, alpha_p, s->x);
        strcpy(s->dist, "post");
        zq += 1.0 / s->d;
        sum_p += s->o;
        sum_q += s->d;
    }
    zq /= (double)(se2 - se1);

    // set the normalization constant for Q(*) points
    for (s = se1; s != se2; ++s)
        s->z = zq;

    zr = zp / zq;
    zri = zq / zp;

    /* Note:
       p(x) := P(x) * zp
       q(x) := Q(x) * zq

       Shrink the exterior region until it does not include any point
       that has P(x) / Q(x) < 1.
       Equivalently:
       (p(x) / zp) / (q(x) / zq) < 1
       p(x) / q(x) < (zp / zq)
     */
    struct ival cut[] = { { 1, 0 }, { 1, 0 }, { 1, 0 }, { 1, 0 } };
    unsigned dim;

    // Find the smallest subregion containing all points in Q(*)
    /* for (s = se1; s != se2; ++s) */
    /*     for (dim = 0; dim != 4; ++dim) */
    /*     { */
    /*         cut[dim].lo = min(cut[dim].lo, s->x[dim]); */
    /*         cut[dim].hi = max(cut[dim].hi, s->x[dim]); */
    /*     } */

    // Expand interior region to include valid points
    unsigned npoints_post_interior = 0;
    sum_p = 0;
    for (s = se1; s != se2; ++s)
        // if (s->o / sum_p > 0.01)
        // if (s->d / s->o < 5.0 * (zq / zp))
        if (s->d / sum_q > 0.001 / (double)npoints_post)
        {
            // include this point
            ++npoints_post_interior;
            sum_p += s->o;
            for (dim = 0; dim != 4; ++dim)
            {
                cut[dim].lo = min(cut[dim].lo, s->x[dim]);
                cut[dim].hi = max(cut[dim].hi, s->x[dim]);
            }
        }

    // Contract interior to exclude too-highly-weighted points
    // This should be done judiciously, just to shrink in one dimension only.
    // Each time a point is excluded in this way, the 
    int more_work = 1;
    unsigned npoints_initial = 1;
    double sum_p_initial;
    while (npoints_initial && npoints_initial != npoints_post_interior)
    {
        npoints_initial = npoints_post_interior;
        npoints_post_interior = 0;
        sum_p_initial = sum_p;
        sum_p = 0;
        for (s = se1; s != se2; ++s)
            if (BETWEEN(s->x[0], cut[0])
                && BETWEEN(s->x[1], cut[1])
                && BETWEEN(s->x[2], cut[2])
                && BETWEEN(s->x[3], cut[3]))
            {
                ++npoints_post_interior;
                sum_p += s->o;
                if (s->o / sum_p_initial > 10.0 / (double)npoints_initial)
                {
                    unsigned imin = min4_index(s->x), imax = max4_index(s->x);
                    if (s->x[imin] < 1.0 - s->x[imax])
                        // hike up low bound to exclude this point
                        cut[imin].lo = max(cut[imin].lo, s->x[imin]);
                    else
                        // move high bound down to exclude this point
                        cut[imax].hi = min(cut[imax].hi, s->x[imax]);
                }
            }
    }
    
    fprintf(stdout, "zp: %g, zq: %g, cut[]: [%g,%g],[%g,%g],[%g,%g],[%g,%g]\n",
            zp, zq,
            cut[0].lo, cut[0].hi, cut[1].lo, cut[1].hi,
            cut[2].lo, cut[2].hi, cut[3].lo, cut[3].hi);

    // initialize pws to point to the union of points:
    // all P(*) points in the periphery and Q(*) points in the
    // interior.
    unsigned i;
    double z_used = 0, z_p = 0, z_q = 0, z_all = 0;
    i = 0;

    unsigned ne = 0, ni = 0; // number of points in exterior or interior violating ratio

    // points from P(*) that reside in exterior
    for (s = wsamples; s != se1; ++s)
    {
        z_p += s->o / s->z;
        if (! (BETWEEN(s->x[0], cut[0])
               && BETWEEN(s->x[1], cut[1])
               && BETWEEN(s->x[2], cut[2])
               && BETWEEN(s->x[3], cut[3])))
        {
            pws[i++] = s;
            z_used += s->o / s->z;
        }
    }

    // points from Q(*) that reside in interior
    for (s = se1; s != se2; ++s)
    {
        z_q += s->o / s->z;
        if (BETWEEN(s->x[0], cut[0])
            && BETWEEN(s->x[1], cut[1])
            && BETWEEN(s->x[2], cut[2])
            && BETWEEN(s->x[3], cut[3]))
        {
            pws[i++] = s;
            z_used += s->o / s->z;
        }
    }

    z_all = z_p + z_q;

    unsigned npoints_used = i;

    // print out complete statistics used to generate CDFs
    fprintf(points_fh, "A\tC\tG\tT\tPointIndex\tratio\tweight\tz-value\tBaseDist\n");
    for (i = 0; i != npoints_used; ++i)
    {
        s = pws[i];
        fprintf(points_fh, "%5.4g\t%5.4g\t%5.4g\t%5.4g\t%u\t%8.4g\t%5.4g\t%5.4g\t%s\n",
                s->x[0], s->x[1], s->x[2], s->x[3], 
                (unsigned)(s - wsamples),
                s->d / s->o, s->o / s->z / z_used, s->z, s->dist);
    }
    
    // print out numerical CDFs
    char base[] = "ACGT";
    char fmt[] = "%s\t%c\t%20.18le\t%u\t%20.18le\t%20.18le\t%g\t%s\n"; 
    char hdr[] = "Dist\tBase\tComp\tPointIndex\tCDF\tPDF\tz-value\tBaseDist\n";
    fprintf(cdf_fh, hdr);

    // exterior points from P(*), interior from Q(*)
    print_cdf(pws, npoints_used, wsamples, "combined", cdf_fh);

    i = 0;
    for (s = wsamples; s != se1; ++s)
        pws[i++] = s;

    // all points in the P(*) sampling
    print_cdf(pws, i, wsamples, "allprior", cdf_fh);

    i = 0;
    for (s = se1; s != se2; ++s)
        pws[i++] = s;

    // all points in the Q(*) sampling
    print_cdf(pws, i, wsamples, "allpost", cdf_fh);

    i = 0;
    for (s = se1; s != se2; ++s)
        ;

    print_cdf(pws, i, wsamples, "intpost", cdf_fh);
        


    // now, for comparison, print out the gold standard
    struct wpoint *gold_points = (struct wpoint *)malloc(npoints_gold * sizeof(struct wpoint));
    struct wpoint **pg = (struct wpoint **)malloc(npoints_gold * sizeof(struct wpoint *));

    
    se1 = gold_points + npoints_gold;
    double alpha_combined[] = { 
        alpha_p[0] + alpha_q[0] - 1,
        alpha_p[1] + alpha_q[1] - 1,
        alpha_p[2] + alpha_q[2] - 1,
        alpha_p[3] + alpha_q[3] - 1
    };

    struct ival region[] = { {0, 1}, {0, 1}, {0, 1}, {0, 1} };
    for (s = gold_points; s != se1; ++s)
    {
        gsl_ran_dirichlet(rand, 4, alpha_combined, s->x);
        // ran_dirichlet_clamped(rand, 4, region, alpha_combined, s->x);
        s->d = gsl_ran_dirichlet_pdf(4, alpha_combined, s->x);
    }

    s = gold_points;
    for (i = 0; s != se1; ++s)
        pg[i++] = s;
    
    for (dim = 0; dim != 4; ++dim)
    {
        cmp_dim = dim;
        qsort(pg, npoints_gold, sizeof(struct wpoints *), pcmp5);
        double sum = 0;
        for (i = 0; i != npoints_gold; ++i)
            fprintf(cdf_fh, fmt, "dgold", base[dim], i, pg[i]->x[dim], 
                    (sum += 1.0 / (double)npoints_gold), pg[i]->d, 1.0, "dgold");
    }


    /* compute component-wise marginal estimates.
     */
    double q[] = { 0.0001, 0.001, 0.01, 0.1, 0.5, 0.9, 0.99, 0.999, 0.9999 };
    double qv[4][sizeof(q) / sizeof(double)];
  
    
    
    // print out marginal estimates
    unsigned iq;
    for (dim = 0; dim != 4; ++dim)
    {
        /* for (iq = 0; iq != sizeof(qv[dim]) / sizeof(double); ++iq) */
        /*     printf(iq ? "\t%f" : "%f", qv[dim][iq]); */
        /* printf("\n"); */
    }

    free(wsamples);
    free(pws);
    free(gold_points);
    free(pg);
    gsl_rng_free(rand);

    fclose(cdf_fh);
    fclose(points_fh);

    return 0;
}
