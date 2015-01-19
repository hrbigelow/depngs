#include "metropolis_sampling.h"
#include <stdio.h>
#include <string.h>
#include <math.h>
#include "yepMath.h"
#include "yepCore.h"

#include <gsl/gsl_statistics_double.h>
#include <gsl/gsl_randist.h>

#define NDIM 4


/*
  Calculates the average of N components of autocorrelation AC[i] =
  Covar(P(i), P(i, offset)) / (SD(P(i)) * SD(P(i, offset))) where P(i)
  is the vector of the i'th component of each point starting at point
  0. P(i, offset) is the vector of the i'th component of each point,
  starting at point 'offset'.
 */
double componentwise_autocorrelation(const double *points, 
                                     size_t num_points, 
                                     size_t offset)
{
    
    double sdev[NDIM], covar[NDIM], mean[NDIM], autocors = 0;
    size_t d;
    for (d = 0; d != NDIM; ++d)
    {
        mean[d] = gsl_stats_mean(points + d, NDIM, num_points);
        sdev[d] = gsl_stats_sd_m(points + d, NDIM, num_points, mean[d]);
        covar[d] = gsl_stats_covariance_m(points + d, NDIM, 
                                          points + (offset * NDIM) + d, NDIM,
                                          num_points - offset,
                                          mean[d], mean[d]);
        autocors += covar[d] / (sdev[d] * sdev[d]);
    }
    autocors /= NDIM;
    return autocors;
}

//returns the lowest offset between samples such that the overall autocorrelation
//at that offset is below <valid_autocor>, or if not found, returns <autocor_max_offset>
size_t best_autocorrelation_offset(const double *samples,
                                   size_t num_samples,
                                   size_t autocor_max_offset,
                                   double valid_autocor)
{
    size_t best_offset = autocor_max_offset;
    double autocor_measure;

    //linear search to find first qualifying offset
    for (size_t o = 1; o != autocor_max_offset && o <= num_samples; ++o)
    {
        autocor_measure = componentwise_autocorrelation(samples, num_samples, o);
        if (autocor_measure < valid_autocor)
        {
            best_offset = o;
            break;
        }
    }

    return best_offset;
}


static inline double alphas_from_counts(const struct packed_counts *cts, 
                                        const double *prior_alpha, 
                                        double *est_alpha)
{
    char base;
    size_t d, qual, strand;
    double alpha0 = prior_alpha[0] + prior_alpha[1]
        + prior_alpha[2] + prior_alpha[3];
        
    memcpy(est_alpha, prior_alpha, sizeof(double) * NDIM);
    for (d = 0; d != cts->num_data; ++d)
    {
        decode_nucleotide(cts->stats_index[d], &base, &qual, &strand);
        alpha0 += cts->stats[d].ct;
        est_alpha[base_to_index(base)] += cts->stats[d].ct;
    }
    return alpha0;
}


static inline void bounded_alphas_from_mean(double *mean, double alpha0, 
                                            const double *bound, double *alphas)
{
    unsigned d;
    double qual_alpha0 = 0, new_alpha0 = 0, adjust;
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

    /* double sum = alphas[0] + alphas[1] + alphas[2] + alphas[3]; */
}



/* */
size_t tune_proposal(const struct packed_counts *cts,
                     const struct posterior_settings *pset,
                     double *proposal_alpha,
                     double *estimated_mean,
                     double *points_buf,
                     struct eval_counts *eval)
{
    double alpha0 = alphas_from_counts(cts, pset->prior_alpha, proposal_alpha);

    size_t i, j, d, cumul_aoff, cur_aoff;

    for (i = 0; i != pset->max_tuning_iterations; ++i)
    {
        cumul_aoff = cur_aoff = 1;
        for (j = 0; j != 3; ++j)
        {
            metropolis_sampling(0, pset->tuning_n_points, cts, pset->logu,
                                proposal_alpha, pset->prior_alpha,
                                cumul_aoff, points_buf, eval);
            
            cur_aoff =
                best_autocorrelation_offset(points_buf,
                                            pset->tuning_n_points,
                                            pset->autocor_max_offset, 
                                            pset->autocor_max_value);
            
            if (cur_aoff == 1) break;
            cumul_aoff *= cur_aoff;
        }
        if (cumul_aoff <= pset->target_autocor_offset) break;
        
        for (d = 0; d != NDIM; ++d)
            estimated_mean[d] = 
                gsl_stats_mean(points_buf + d, NDIM, pset->tuning_n_points);

        bounded_alphas_from_mean(estimated_mean, alpha0, 
                                 pset->prior_alpha, proposal_alpha);
    }
    return cumul_aoff;
}

/* Generate sample points according to Metropolis criterion, using the
    alphas for a dirichlet proposal distribution, and the polynomial
    terms for the posterior. The sample points generated are written
    to sample_points, in the range [start_point, end_point).  Collect
    every 'nth' point.  Note the following identities, and the
    cancellation of the prior.

   Post(pri,obs) := Dir(pri)Mult(obs)
   Prop(pro) := Dir(pri+(pro - pri)) := Dir(pri)Dir(pro-pri+1)
   Post(pri,obs) / Prop(pri+obs) = Mult(obs) / Dir(pro-pri+1)

   let res = pro-pri+1  (residual alphas)
*/
#define BATCH 16

void
metropolis_sampling(unsigned short start_point, 
                    unsigned short end_point,
                    const struct packed_counts *cts,
                    const double *logu, /* must have start_point
                                           + end_point points */
                    double *proposal_alpha,
                    const double *prior_alpha,
                    size_t nth, /* collect every nth point */
                    double *sample_points,
                    struct eval_counts *eval)
{
    const struct cpd_count *trm = cts->stats, *trm_end = trm + cts->num_data;

    double *pa = proposal_alpha;
    fprintf(stdout, "metropolis_sampling: %f,%f,%f,%f\n", pa[0], pa[1], pa[2], pa[3]);
    fflush(stdout);

    double points[BATCH * NDIM], dotp[BATCH], ldotp[BATCH];
    double llh[BATCH], rll[BATCH], ivp[BATCH];
    double ar; /* acceptance ratio */
    double *pt, *outpt = sample_points + start_point * NDIM;
    size_t i;
    gsl_rng *randgen = gsl_rng_alloc(gsl_rng_taus);

    unsigned short a = start_point, u = start_point;

    memset(llh, 0, sizeof(llh));

    double residual_alpha[] = {
        proposal_alpha[0] - prior_alpha[0] + 1,
        proposal_alpha[1] - prior_alpha[1] + 1,
        proposal_alpha[2] - prior_alpha[2] + 1,
        proposal_alpha[3] - prior_alpha[3] + 1
    };

    /* there are end_point values in logu.  and, in one iteration of
     this while-loop, we will increment u BATCH * nth times. So, this
     is the max value u can have before the iteration. */
    unsigned u_max = end_point - nth * BATCH;
    while (a != end_point)
    {
        if (u >= u_max) u = start_point;

        /* 1. Generate a batch of proposal points and log likelihoods
           for them */
        pt = points;
        for (i = 0; i != BATCH; ++i)
        {
            gsl_ran_dirichlet(randgen, NDIM, proposal_alpha, pt);
            rll[i] = gsl_ran_dirichlet_lnpdf(NDIM, residual_alpha, pt);
            pt += NDIM;
        }        
        eval->n_dirichlet += BATCH;

        /* 3. Generate a batch log likelihoods */
        while (trm != trm_end)
        {
            (void)yepCore_DotProduct_V64fV64f_S64f(points, trm->cpd, dotp, NDIM);
            (void)yepMath_Log_V64f_V64f(dotp, ldotp, BATCH);
            (void)yepCore_Multiply_IV64fS64f_IV64f(ldotp, (double)trm->ct, BATCH);
            (void)yepCore_Add_V64fV64f_V64f(llh, ldotp, llh, BATCH);
            ++trm;
        }
        eval->n_yep += cts->num_data * BATCH;

        /* 4. Compute ivp ratios */
        for (i = 0; i != BATCH; ++i)
            ivp[i] = rll[i] - llh[i]; /* any way to avoid the division? */
        
        /* 5. Copy with rejection */
        unsigned short p = 0, c;
        for (c = 1; c != BATCH && a != end_point; ++c)
        {
            ar = ivp[c] - ivp[p];
            if (ar > logu[++u]) p = c;
            if (u % nth == 0)
            {
                memcpy(outpt, &points[p * NDIM], NDIM * sizeof(double));
                ++a;
                outpt += NDIM;
            }
        }
    }
    gsl_rng_free(randgen);
}


