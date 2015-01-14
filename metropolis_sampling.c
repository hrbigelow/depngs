#include "metropolis_sampling.h"
#include <stdio.h>
#include <string.h>
#include <math.h>
#include "yepMath.h"

#define BATCH 16
#define NDIM 4


inline double alphas_from_counts(struct packed_counts *cts, 
                                 double *prior_alpha, double *est_alpha)
{
    char base;
    size_t d, qual, strand;
    double alpha0 = 0;
    memcpy(est_alpha, prior_alpha, sizeof(double) * NDIM);
    for (d = 0; d != cts->num_data; ++d)
    {
        decode_nucleotide(cts->stats_index[d], &base, &qual, &strand);
        alpha0 += cts->stats[d].ct;
        est_alpha[base_to_index(base)] += cts->stats[d].ct;
    }
    return alpha0;
}


inline void bounded_alphas_from_mean(double *mean, double alpha0, 
                                     double *bound, double *alphas)
{
    unsigned d;
    double qual_alpha0 = 0, new_alpha0 = 0, sum, adjust;
    for (d = 0; d != NDIM; ++d)
    {
        alphas[d] = mean[d] * alpha0;
        if (alphas[d] < bound[d]) alphas[d] = bound[d];
        else qual_alpha0 += alphas[d];
        new_alpha0 += alphas[d];
    }

    adjust = 1.0 + (alpha0 - new_alpha0) / qual_alpha0;
    for (d = 0; d != NDIM; ++d)
        if (alphas[d] != bound[d]) alphas[d] *= adjust;

    sum = alphas[0] + alphas[1] + alphas[2] + alphas[3];
}



/* */
size_t tune_proposal(struct packed_counts *cts,
                     const struct posterior_settings *set,
                     double *proposal_alpha,
                     double *estimated_mean,
                     double *points_buf)
{
    double alpha0 = alphas_from_counts(cts, set->prior_alpha, proposal_alpha);

    size_t i, j, cumul_aoff, cur_aoff;
    for (i = 0; i != set->max_tuning_iterations; ++i)
    {
        cumul_aoff = current_aoff = 1;
        for (j = 0; j != 3; ++j)
        {
            metropolis_sampling(0, set->tuning_n_points, cts->stats, cts->num_data,
                                proposal_alpha, logu, cumul_aoff, points_buf);
            
            current_aoff =
                best_autocorrelation_offset(sample_points_buf, NDIM,
                                            set->tuning_n_points,
                                            set->autocor_max_offset, 
                                            set->autocor_max_value);
            
            if (current_aoff == 1) break;
            cumul_aoff *= current_aoff;
        }
        if (cumul_aoff <= set->target_autocor_offset) break;
        
        multivariate_mean(points_buf, NDIM, set->tuning_n_points, estimated_mean);
        bounded_alphas_from_mean(estimated_mean, alpha0, prior_alpha, proposal_alpha);
    }
    return cumul_aoff;
}


/* Generate sample points according to Metropolis criterion, using the
    alphas for a dirichlet proposal distribution, and the polynomial
    terms for the posterior. The sample points generated are written
    to sample_points, in the range from start_point to
    n_points_wanted.  Collect every 'nth' point */
void metropolis_sampling(unsigned short start_point, unsigned short n_points_wanted,
                         const struct packed_counts *cts,
                         const double *logu, /* logs of U[0, 1] values */
                         double *proposal_alpha,
                         size_t nth, /* collect every nth point */
                         double *sample_points)
{
    const struct cpd_count *trm_end = cts->stats + cts->num_data;

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
            gsl_ran_dirichlet(proposal_alpha, NDIM, pt);
            rll[i] = gsl_ran_dirichlet_lnpdf(proposal_alpha, NDIM, pt);
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


