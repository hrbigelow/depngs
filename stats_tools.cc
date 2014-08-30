#include "stats_tools.h"

#include <cstdio>
#include <numeric>
#include <functional>
#include <cmath>
#include <algorithm>


#include <gsl/gsl_statistics_double.h>
#include <gsl/gsl_math.h>

void multivariate_mean(const double *samples,
                       size_t ndim,
                       size_t num_samples,
                       double *mean)
{
    //update the mean's
    for (size_t d = 0; d != ndim; ++d)
        mean[d] = gsl_stats_mean(samples + d, ndim, num_samples);
}


void multivariate_mean_covariance(const double *samples,
                                  size_t ndim,
                                  size_t num_samples,
                                  double *mean, double *covariance)
{

    //update the mean's
    for (size_t d = 0; d != ndim; ++d)
        mean[d] = gsl_stats_mean(samples + d, ndim, num_samples);
    
    //Update the covariances
    for (size_t d1 = 0; d1 != ndim; ++d1)
    {
        for (size_t d2 = 0; d2 != ndim; ++d2)
        {
            covariance[d1 * ndim + d2] = 
                gsl_stats_covariance_m(samples + d1, ndim,
                                       samples + d2, ndim,
                                       num_samples, mean[d1], mean[d2]);
        }
    }
}


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
    
    double sdev[NUM_NUCS], covar[NUM_NUCS], mean[NUM_NUCS];
    double autocors = 0;

    multivariate_mean(points, NUM_NUCS, num_points, mean);

    for (size_t d = 0; d != NUM_NUCS; ++d)
    {
        sdev[d] = gsl_stats_sd_m(points + d, NUM_NUCS, num_points, mean[d]);
        covar[d] = gsl_stats_covariance_m(points + d, NUM_NUCS, 
                                          points + (offset * NUM_NUCS) + d, NUM_NUCS,
                                          num_points - offset,
                                          mean[d], mean[d]);
        autocors += covar[d] / (sdev[d] * sdev[d]);
    }
    autocors /= NUM_NUCS;
    return autocors;
}



//returns the lowest offset between samples such that the overall autocorrelation
//at that offset is below <valid_autocor>, or if not found, returns <autocor_max_offset>
size_t best_autocorrelation_offset(const double *samples,
                                   size_t ndim,
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


bool all_positive(const double *x, size_t n)
{
    for (size_t d = 0; d != n; ++d)
        if (x[d] < 0)
            return false;

    return true;
}

bool normalized(const double *x, size_t n, double delta)
{
    double sum = 0.0;
    for (size_t d = 0; d != n; ++d)
        sum += x[d];

    return gsl_fcmp(sum, 1.0, delta) == 0;
}


void normalize(double *x, size_t n, double *x_out)
{
    double sum = std::accumulate(x, x + n, 0.0);
    std::transform(x, x + n, x_out, std::bind2nd(std::divides<double>(), sum));
}



//calculate D_kl(P||Q) see
//http://en.wikipedia.org/wiki/Kullback-Leibler_divergence 

//assume that q[n] > 0 || q[n] == 0 && p[n] == 0 this is consistent
//with the use of q as the frequency of a subsampling from p
double kullback_leibler_divergence(const double *p, const double *q, size_t N)
{
    double kl = 0.0;

    for (size_t n = 0; n != N; ++n)
    {
        if (p[n] > 0.0)
        {
            if (q[n] == 0.0)
            {
                fprintf(stderr, "Error: KL divergence assumes "
                        "q[n] > 0 || q[n] == 0 && p[n] == 0\n");
                exit(1);
            }
            else
            {
                kl += p[n] * gsl_sf_log(p[n]/q[n]);
            }
        }
        else
        {
        }
    }
    return kl;
}
