#include "dirichlet_points_gen.h"
#include "nucleotide_stats.h"

#include <gsl/gsl_randist.h>
#include <gsl/gsl_sf_exp.h>

#include <string.h>
#include "yepMath.h"
#include "yepCore.h"


/* compute GEN_POINTS_BATCH of unnormalized dirichlet values using the
   same alpha.  assume NUM_NUCS dimensions */
void ran_dirichlet_lnpdf_unnormalized(double *alpha, double *points, double *lndir)
{
    double 
        lnpoints[GEN_POINTS_BATCH * NUM_NUCS], 
        dir_lnpdf[GEN_POINTS_BATCH * NUM_NUCS],
        alpha_minus_one[NUM_NUCS], 
        *lnp;

    size_t i, p;
    for (i = 0; i != NUM_NUCS; ++i) alpha_minus_one[i] = alpha[i] - 1.0;

    (void)yepMath_Log_V64f_V64f(points, lnpoints, sizeof(lnpoints) / sizeof(lnpoints[0]));

    lnp = lnpoints;
    for (p = 0; p != GEN_POINTS_BATCH; ++p)
    {
        dir_lnpdf[p] = 0.0;
        for (i = 0; i != NUM_NUCS; ++i, ++lnp)
            dir_lnpdf[p] += alpha_minus_one[i] * *lnp;
    }
    memcpy(lndir, dir_lnpdf, sizeof(dir_lnpdf));
}


#if 0
void gsl_ran_dirichlet_batched(const gsl_rng *r, 
                               const double *alpha, double *theta)
{
    size_t i;
    double norm = 0.0, norm_inv;

    for (i = 0; i != NUM_NUCS; i++)
        theta[i] = gsl_ran_gamma(r, alpha[i], 1.0);
    
    for (i = 0; i != NUM_NUCS; i++)
        norm += theta[i];
    
    if (norm < GSL_SQRT_DBL_MIN)  /* Handle underflow */
    {
        ran_dirichlet_small(r, NUM_NUCS, alpha, theta);
        return;
    }

    norm_inv = 1.0 / norm;
    for (i = 0; i != NUM_NUCS; i++) theta[i] *= norm_inv;
}
#endif



/* the following four functions can be used with binomial_est's struct
   points_gen. */
void gen_dirichlet_points_wrapper(const void *par, POINT *points)
{
    int i;
    struct dir_points_par *gd = (struct dir_points_par *)par;

    for (i = 0; i != GEN_POINTS_BATCH; ++i)
    {
        gsl_ran_dirichlet(gd->randgen, NUM_NUCS, gd->alpha, *points);
        ++points;
    }
}


/* generate a 'reference' point, representing the corner of the
   simplex corresponding to the reference base, or a point outside the
   simplex for reference 'N'.  This external point will serve as an
   'always different' point. */
void gen_reference_points_wrapper(const void *par, POINT *points)
{
    static POINT ref_points[] = {
        { 1, 0, 0, 0 },
        { 0, 1, 0, 0 },
        { 0, 0, 1, 0 },
        { 0, 0, 0, 1 },
        { -1, -1, -1, -1 }
    };
    char *refbase = (char *)par;
    int ref_ind = base_to_index(*refbase);
    
    int i;
    for (i = 0; i != GEN_POINTS_BATCH; ++i, ++points)
        memcpy(*points, ref_points[ref_ind], sizeof(ref_points[0]));
}


void calc_post_to_dir_ratio(POINT *points, const void *par, double *weights)
{
    int i;
    struct calc_post_to_dir_par *pd = (struct calc_post_to_dir_par *)par;
    const struct cpd_count 
        *trm = pd->post_counts->stats, 
        *trm_end = trm + pd->post_counts->num_data;

    double 
        dotp[GEN_POINTS_BATCH], ldotp[GEN_POINTS_BATCH],
        llh[GEN_POINTS_BATCH], ldir[GEN_POINTS_BATCH];

    POINT *p, *pe = points + GEN_POINTS_BATCH;
    while (trm != trm_end)
    {
        for (p = points, i = 0; p != pe; ++p, ++i)
            dotp[i] = 
                (*p)[0] * trm->cpd[0] + (*p)[1] * trm->cpd[1]
                + (*p)[2] * trm->cpd[2] + (*p)[3] * trm->cpd[3];
        
        (void)yepMath_Log_V64f_V64f(dotp, ldotp, GEN_POINTS_BATCH);
        (void)yepCore_Multiply_IV64fS64f_IV64f(ldotp, (double)trm->ct, GEN_POINTS_BATCH);
        (void)yepCore_Add_V64fV64f_V64f(llh, ldotp, llh, GEN_POINTS_BATCH);
        ++trm;
    }
    ran_dirichlet_lnpdf_unnormalized(pd->proposal_alpha, (double *)points, ldir);

    for (i = 0; i != GEN_POINTS_BATCH; ++i)
        weights[i] = gsl_sf_exp(llh[i] - ldir[i]);
}

void calc_dummy_ratio(POINT *point, const void *par, double *weights)
{
    int i;
    for (i = 0; i != GEN_POINTS_BATCH; ++i)
        weights[i] = 1;
}



