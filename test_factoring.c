#include <gsl/gsl_randist.h>
#include <gsl/gsl_cdf.h>
#include <gsl/gsl_sf_log.h>
#include <gsl/gsl_sf_exp.h>
#include <stdio.h>
#include <float.h>
#include <string.h>
#include <math.h>
// #include <algorithm>
#include <assert.h>

struct wpoint {
    double x[4], d, o, ln_d, ln_o;
    double ln_cumul_o; // density and 'other'.
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


// compare wpoint *pa and *pb
int ocmp(const void *pa, const void *pb)
{
    const struct wpoint 
        *a = *(struct wpoint **)pa, 
        *b = *(struct wpoint **)pb;

    return a->ln_o < b->ln_o ? -1 : (a->ln_o > b->ln_o ? 1 : 0);
}


void *lower_bound(void *pbase, size_t memb, size_t size,
                  int (comp)(const void *, const void *),
                  void *val)
{
    char *mid, *base = (char *)pbase;
    size_t half;
    while (memb > 0)
    {
        half = memb >> 1;
        mid = base + (size * half);
        if (comp((void *)mid, val) == -1)
        {
            base = mid;
            base += size;
            memb -= half + 1;
        }
        else
            memb = half;
    }
    return (void *)base;
}


inline double max(double x, double y) { return x > y ? x : y; }
inline double min(double x, double y) { return x < y ? x : y; }

#define MIN4(a) min(min((a)[0], (a)[1]), min((a)[2], (a)[3]))
#define MAX4(a) max(max((a)[0], (a)[1]), max((a)[2], (a)[3]))


/* compute ln(a + b) safely, given ln(a) and ln(b). */
double safe_sum(double ln_a, double ln_b)
{
    if (isnormal(ln_a) && isnormal(ln_b))
        return ln_a < ln_b 
            ? ln_b + gsl_sf_log(gsl_sf_exp(ln_a - ln_b) + 1.0)
            : ln_a + gsl_sf_log(gsl_sf_exp(ln_b - ln_a) + 1.0);
    else if (isnormal(ln_a))
        return ln_a;
    else 
        return ln_b;
}

// print out the CDF given a number of points.
// normalize to the given constant
void print_cdf(struct wpoint **points, unsigned npoints,
               double total_mass,
               struct wpoint *first_point, // needed to give each point an index
               const char *dist_name,
               FILE *cdf_fh)
{
    unsigned dim;
    char fmt[] = "%s\t%c\t%20.18le\t%u\t%20.18le\t%g\t%s\n"; 
    char base[] = "ACGT";
    struct wpoint *s, **sp, **se = points + npoints;

    double ln_z = gsl_sf_log(total_mass);

    double ln_sum, max_ln_step = -INFINITY;
    for (dim = 0; dim != 4; ++dim)
    {
        cmp_dim = dim;
        qsort(points, npoints, sizeof(struct wpoints *), pcmp5);
        ln_sum = NAN;

        for (sp = points; sp != se; ++sp)
        {
            s = *sp;
            max_ln_step = max(max_ln_step, s->ln_o - ln_z);
            ln_sum = safe_sum(ln_sum, s->ln_o);
                fprintf(cdf_fh, 
                        fmt, dist_name, base[dim], (unsigned)(s - first_point),
                        s->x[dim], gsl_sf_exp(ln_sum - ln_z), s->d, s->dist);
        }
    }
    fprintf(stderr, "%s, max_ln_step = %g, max_step = %g\n", dist_name, max_ln_step, gsl_sf_exp(max_ln_step));

}


void populate_points(struct wpoint *points, unsigned npoints,
                     double *alpha_dist, double *alpha_weights,
                     gsl_rng *rand,
                     const char *name,
                     double *dist_norm)
{
    struct wpoint *s, *se = points + npoints;
    double ln_dist_norm = NAN;
    double dn = 0;
    for (s = points; s != se; ++s)
    {
        gsl_ran_dirichlet(rand, 4, alpha_dist, s->x);
        s->d = gsl_ran_dirichlet_pdf(4, alpha_dist, s->x);
        s->o = gsl_ran_dirichlet_pdf(4, alpha_weights, s->x);
        s->ln_d = gsl_ran_dirichlet_lnpdf(4, alpha_dist, s->x);
        s->ln_o = gsl_ran_dirichlet_lnpdf(4, alpha_weights, s->x);
        strcpy(s->dist, name);
        dn += 1.0 / s->d;
        ln_dist_norm = safe_sum(ln_dist_norm, -s->ln_d);
    }
    ln_dist_norm -= gsl_sf_log(npoints);
    dn /= (double)npoints;
    *dist_norm = gsl_sf_exp(ln_dist_norm);
}


/* normalize densities.  after this, sum(1/d) should be 1, and sum(o)
   should be the same for any distribution */
void normalize_points(struct wpoint *points, unsigned npoints,
                      double d_norm, double o_norm)
{
    struct wpoint *s, *se = points + npoints;
    double ln_d_norm = gsl_sf_log(d_norm);
    double ln_o_norm =  gsl_sf_log(o_norm);
    for (s = points; s != se; ++s)
    {
        s->ln_d -= ln_d_norm;
        s->ln_o -= ln_o_norm;
        s->d /= d_norm;
        s->o /= o_norm;
    }
}


double sum_of_weights(struct wpoint **pws, unsigned npoints)
{
    struct wpoint **p, **pe = pws + npoints;
    double sum = 0;
    for (p = pws; p != pe; ++p)
        sum += (*p)->o;
    return sum;
}

/* find ln_thresh such that all points with ln_o < ln_thresh are
   'smooth' */
double smooth_threshold(struct wpoint **points, double mult, unsigned *npoints)
{
    unsigned *n = npoints;
    struct wpoint **pp, **pe = points + *n;
    double ln_total = NAN;
    for (pp = points; pp != pe; ++pp)
        ln_total = safe_sum(ln_total, (*pp)->ln_o);

    qsort(points, *n, sizeof(struct wpoint *), ocmp);
    struct wpoint query, *pq = &query, **lb;
    
    /* the log form of mult * avg_over_n(o) */
    query.ln_o = gsl_sf_log(mult) + ln_total - gsl_sf_log(*n);

    lb = lower_bound(points, *n, sizeof(struct wpoint *), ocmp, &pq);
    *n = lb - points;
    
    return lb == pe ? DBL_MAX : (*lb)->ln_o;
}


double max_to_avg_ratio(struct wpoint **points, unsigned npoints)
{
    struct wpoint **pp, **pe = points + npoints;
    double ln_total = NAN, maxval = -DBL_MAX;
    for (pp = points; pp != pe; ++pp)
    {
        ln_total = safe_sum(ln_total, (*pp)->ln_o);
        maxval = max(maxval, (*pp)->o);
    }
    return maxval / (gsl_sf_exp(ln_total) / (double)npoints);
}

int main(int argc, char **argv)
{
    char *str_alpha_p = argv[1], *str_alpha_q = argv[2];
    char *cdf = argv[3], *points = argv[4];

    unsigned nprior = atoi(argv[5]);
    unsigned npost = atoi(argv[6]);
    unsigned ngold = atoi(argv[7]);
    double step_mult = atof(argv[8]);

    gsl_set_error_handler_off();

    FILE *cdf_fh = fopen(cdf, "w");
    FILE *points_fh = fopen(points, "w");

    double alpha_p[4], alpha_q[4];
    sscanf(str_alpha_p, "%lf,%lf,%lf,%lf", &alpha_p[0], &alpha_p[1], &alpha_p[2], &alpha_p[3]);
    sscanf(str_alpha_q, "%lf,%lf,%lf,%lf", &alpha_q[0], &alpha_q[1], &alpha_q[2], &alpha_q[3]);

    double alpha_combined[] = { 
        alpha_p[0] + alpha_q[0] - 1,
        alpha_p[1] + alpha_q[1] - 1,
        alpha_p[2] + alpha_q[2] - 1,
        alpha_p[3] + alpha_q[3] - 1
    };

    unsigned ntotal = nprior + npost + ngold;
    struct wpoint *buf = (struct wpoint *)malloc(ntotal * sizeof(struct wpoint));
    struct wpoint *p_points = buf, *q_points = buf + nprior, *g_points = q_points + npost;
    struct wpoint **tmp_points = (struct wpoint **)malloc((nprior + npost) * sizeof(struct wpoint *));
    struct wpoint **int_smooth_points = (struct wpoint **)malloc(npost * sizeof(struct wpoint *));
    struct wpoint **ext_points = (struct wpoint **)malloc(nprior * sizeof(struct wpoint *));
    struct wpoint **ext_smooth_points = (struct wpoint **)malloc(nprior * sizeof(struct wpoint *));
    struct wpoint **ext_rough_points = (struct wpoint **)malloc(nprior * sizeof(struct wpoint *));
    struct wpoint **all_smooth_points = (struct wpoint **)malloc((nprior + npost) * sizeof(struct wpoint *));

    struct wpoint **pg = (struct wpoint **)malloc(ngold * sizeof(struct wpoint *));

    struct wpoint **pp, **pp2, **pe;
    struct wpoint *s, *p_end = p_points + nprior, *q_end = q_points + npost;

    gsl_rng *rand = gsl_rng_alloc(gsl_rng_taus);

    // populate P(*)
    double p_norm;
    populate_points(p_points, nprior, alpha_p, alpha_q, rand, "prior", &p_norm);

    // populate Q(*)
    double q_norm;
    populate_points(q_points, npost, alpha_q, alpha_p, rand, "post", &q_norm);

    /* populate gold(*) */
    double g_norm;
    double alpha_uniform[] = { 1, 1, 1, 1 };
    populate_points(g_points, ngold, alpha_combined, alpha_uniform, rand, "gold", &g_norm);

    /* normalize */
    normalize_points(p_points, nprior, p_norm, q_norm);
    normalize_points(q_points, npost, q_norm, p_norm);
    normalize_points(g_points, ngold, g_norm, g_norm);

    unsigned n;
    
    /* gather all Q(*) points (use int_smooth_points temporarily) */
    for (s = q_points, pp = int_smooth_points; s != q_end; ++s)
        *pp++ = s;
    n = pp - int_smooth_points;

    /* collect Q(*) smooth (interior) points */
    double ln_p_bound = smooth_threshold(int_smooth_points, step_mult, &n);
    for (s = q_points, pp = int_smooth_points; s != q_end; ++s)
        if (s->ln_o <= ln_p_bound)
            *pp++ = s;
    unsigned n_int_smooth_q = pp - int_smooth_points;

    double ma_ratio = max_to_avg_ratio(int_smooth_points,);
    double int_mass_on_q = sum_of_weights(int_smooth_points, n_int_smooth_q);


    /* collect exterior P(*) points */
    for (s = p_points, pp = ext_points; s != p_end; ++s)
        if (s->ln_d > ln_p_bound)
            *pp++ = s;
    unsigned n_ext_p = n = pp - ext_points;
    double ln_q_bound = smooth_threshold(ext_points, step_mult, &n);

    /* collect exterior P(*) smooth points */
    for (s = p_points, pp = ext_smooth_points; s != p_end; ++s)
        if (s->ln_d > ln_p_bound && s->ln_o < ln_q_bound)
            *pp++ = s;
    unsigned n_ext_smooth_p = pp - ext_smooth_points;

    
    double ma2_ratio = max_to_avg_ratio(ext_smooth_points, n_ext_smooth_p);

    double ext_mass_on_p = sum_of_weights(ext_smooth_points, n_ext_smooth_p);

    double total_mixed_mass = int_mass_on_q + ext_mass_on_p;

    fprintf(stderr, "Number of points unaccounted for: %u\n", n_ext_p - n_ext_smooth_p);

    // print out numerical CDFs
    char hdr[] = "Dist\tBase\tComp\tPointIndex\tCDF\tPDF\tz-value\tBaseDist\n";
    fprintf(cdf_fh, hdr);

    print_cdf(int_smooth_points, n_int_smooth_q, total_mixed_mass, buf, "int_smooth_q", cdf_fh);
    print_cdf(ext_smooth_points, n_ext_smooth_p, total_mixed_mass, buf, "ext_smooth_p", cdf_fh);

   
    pe = int_smooth_points + n_int_smooth_q;
    for (pp = int_smooth_points, pp2 = all_smooth_points; pp != pe; ++pp, ++pp2)
        *pp2 = *pp;

    pe = ext_smooth_points + n_ext_smooth_p;
    for (pp = ext_smooth_points; pp != pe; ++pp, ++pp2)
        *pp2 = *pp;

    unsigned n_smooth = n_int_smooth_q + n_ext_smooth_p;
    print_cdf(all_smooth_points, n_smooth, total_mixed_mass, p_points, "smooth_combined", cdf_fh);

    /* all points in the P(*) sampling */
    for (s = p_points, pp = tmp_points; s != p_end; ++s)
        *pp++ = s;
    print_cdf(tmp_points, nprior, total_mixed_mass, p_points, "all_p", cdf_fh);

    /* all points in the Q(*) sampling */
    for (s = q_points, pp = tmp_points; s != q_end; ++s)
        *pp++ = s;

    print_cdf(tmp_points, npost, total_mixed_mass, q_points, "all_q", cdf_fh);


    // print out complete statistics used to generate CDFs
    /* fprintf(points_fh, "A\tC\tG\tT\tPointIndex\tratio\tweight\tz-value\tBaseDist\n"); */
    /* for (i = 0; i != npoints_used; ++i) */
    /* { */
    /*     s = pws[i]; */
    /*     fprintf(points_fh, "%5.4g\t%5.4g\t%5.4g\t%5.4g\t%u\t%8.4g\t%5.4g\t%5.4g\t%s\n", */
    /*             s->x[0], s->x[1], s->x[2], s->x[3],  */
    /*             (unsigned)(s - wsamples), */
    /*             s->d / s->o, s->o / s->z / z_used, s->z, s->dist); */
    /* } */
    

    /* print out the gold standard */
    struct wpoint *gold_end = g_points + ngold;
    for (s = g_points, pp = pg; s != gold_end; ++s)
        *pp++ = s;
    
    double gold_mass = sum_of_weights(pg, ngold);
    print_cdf(pg, ngold, gold_mass, g_points, "dgold", cdf_fh);

    /* compute component-wise marginal estimates. */
    /* double q[] = { 0.0001, 0.001, 0.01, 0.1, 0.5, 0.9, 0.99, 0.999, 0.9999 }; */
    /* double qv[4][sizeof(q) / sizeof(double)]; */
    
    // print out marginal estimates
    /* unsigned iq; */
    /* for (dim = 0; dim != 4; ++dim) */
    /* { */
        /* for (iq = 0; iq != sizeof(qv[dim]) / sizeof(double); ++iq) */
        /*     printf(iq ? "\t%f" : "%f", qv[dim][iq]); */
        /* printf("\n"); */
    /* } */

    free(buf);
    free(tmp_points);
    free(int_smooth_points);
    free(ext_points);
    free(ext_smooth_points);
    free(ext_rough_points);
    free(all_smooth_points);
    free(pg);
    gsl_rng_free(rand);

    fclose(cdf_fh);
    fclose(points_fh);

    return 0;
}
