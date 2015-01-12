#include "metropolis_sampling.h"
#include <stdio.h>
#include <string.h>
#include <math.h>


#define BATCH 16;
#define NDIM 4;

/* Generate sample points according to Metropolis criterion, using the
    alphas for a dirichlet proposal distribution, and the polynomial
    terms for the posterior. The sample points generated are written
    to sample_points, in the range from start_point to
    n_points_wanted.  Collect every 'nth' point */
void metropolis_sampling(unsigned short start_point, unsigned short n_points_wanted,
                         const struct cpd_count *term, size_t n_term,
                         double *proposal_alpha,
                         double *prior_alpha,
                         double *logu, /* logs of U[0, 1] values */
                         size_t nth, /* collect every nth point */
                         double *sample_points)
{
    enum YepStatus status = yepLibrary_Init();
    assert(status == YepStatusOk);

    double alpha_residual[] = { 
        proposal_alpha[0] - prior_alpha[0] + 1, 
        proposal_alpha[1] - prior_alpha[1] + 1,
        proposal_alpha[2] - prior_alpha[2] + 1,
        proposal_alpha[3] - prior_alpha[3] + 1 
    };
    
    const struct cpd_count *trm_end = trm + n_trm;

    double points[BATCH * NDIM], dotp[BATCH], ldotp[BATCH];
    double llh[BATCH], rll[BATCH], ivp[BATCH];
    double ar; /* acceptance ratio */
    double *pt, *outpt = sample_points + start_point * NDIM;
    unsigned short a = start_point;

    while (a != n_points_wanted)
    {
        /* 1. Generate a batch of proposal points and log likelihoods for them */
        pt = points;
        for (i = 0; i != BATCH; ++i)
        {
            gsl_ran_dirichlet(alpha_residual, NDIM, pt);
            rll[i] = gsl_ran_dirichlet_lnpdf(alpha_residual, NDIM, pt);
            pt += NDIM;
        }        
        /* 3. Generate a batch log likelihoods */
        while (trm != trm_end)
        {
            (void)yepCore_DotProduct_V64fV64f_S64f(points, trm->cpd, dotp, NDIM);
            (void)yepMath_Log_V64f_V64f(dotp, ldotp, BATCH);
            (void)yepCore_Multiply_IV64fS64f_IV64f(ldotp, (double)trm->ct, BATCH);
            (void)yepCore_Add_V64fV64f_V64f(llh, ldotp, llh, BATCH);
            ++trm;
        }
        /* 4. Compute ivp ratios */
        for (i = 0; i != BATCH; ++i)
            ivp[i] = rll[i] / llh[i]; /* any way to avoid the division? */
        
        /* 5. Copy with rejection */
        unsigned short c, p = 0, u = 0;
        for (c = 1; c != BATCH && a != n_points_wanted; ++c)
        {
            ar = ivp[c] / ivp[p];
            if (ar > logU[++u]) p = c;
            if (u % nth)
            {
                memcpy(outpt, &points[p * NDIM], NDIM * sizeof(double));
                ++a;
                outpt += NDIM;
            }
        }
    }
}
