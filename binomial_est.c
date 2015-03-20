#include "binomial_est.h"

#include <gsl/gsl_cdf.h>
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_errno.h>

#include <stdlib.h>
#include <assert.h>
#include <stdio.h>

#include "cache.h"

/*
Consider generating the distance distribution between a real
locus (a posterior in the simplex) and a 'fake' locus that
represents homozygous reference (a point-mass at one of the
simplex corner points).

Secondly, consider the 1D distribution which is that same real
locus posterior marginalized to that corner point. This marginal
distribution is the same thing as the distance distribution.

So, an easy way to write the function is simply to calculate the
distance distribution.  (This is a bit more expensive, but
perhaps worth it in the benefits of simplifying the code.

*/

#define CDF_ERROR(reason, gsl_errno) GSL_ERROR_VAL(reason, gsl_errno, GSL_NAN)

static double bisect(double x, double P, double a, double b, double xtol, double Ptol)
{
    double x0 = 0, x1 = 1, Px;
    while (fabs(x1 - x0) > xtol) {
        Px = gsl_cdf_beta_P (x, a, b);
        if (fabs(Px - P) < Ptol) return x;  
        else if (Px < P) x0 = x;
        else if (Px > P) x1 = x;
        x = 0.5 * (x0 + x1);
    }
    return x;
}  

double beta_Qinv(double P, double a, double b);

/* hack from GSL that doesn't give 'fail to converge' error.
   same as gsl_cdf_beta*/
double beta_Pinv(double P, double a, double b)
{
    double x, mean;

    if (P < 0.0 || P > 1.0) CDF_ERROR("P must be in range 0 < P < 1", GSL_EDOM);

    if (a < 0.0) CDF_ERROR("a < 0", GSL_EDOM);

    if (b < 0.0) CDF_ERROR("b < 0", GSL_EDOM);

    if (P == 0.0) return 0.0;

    if (P == 1.0) return 1.0;

    if (P > 0.5) return beta_Qinv(1 - P, a, b);

    mean = a / (a + b);

    if (P < 0.1)
    {
        /* small x */
        double lg_ab = gsl_sf_lngamma(a + b);
        double lg_a = gsl_sf_lngamma(a);
        double lg_b = gsl_sf_lngamma(b);
      
        double lx = (log(a) + lg_a + lg_b - lg_ab + log(P)) / a;
        if (lx <= 0) {
            x = exp(lx);                    /* first approximation */
            x *= pow(1 - x, -(b - 1) / a);  /* second approximation */
        } else
            x = mean;
      
        if (x > mean) x = mean;
    }
    else x = mean; /* Use expected value as first guess */
  
    /* Do bisection to get to within tolerance */
    x = bisect(x, P, a, b, 0.01, 1e-10);
    return x;

}


double beta_Qinv(double Q, double a, double b)
{
    if (Q > 0.5) return beta_Pinv(1 - Q, a, b);
    else return 1 - beta_Pinv(Q, b, a);
}


/* Maximum N of precalculated Beta values. */
#define NUM_BETA_PRECALC 1000

/* beta_lo[f][n] = gsl_cdf_beta_Pinv(beta_conf, s + 1/2, n - s + 1/2), where
   s: # of successes
   f: # of failures (n - s)
   n: # of trials
   This is the low bound of the Jeffrey's posterior for binomial
   confidence intervals. */
double beta_lo[NUM_BETA_PRECALC][NUM_BETA_PRECALC];


/* beta_hi[f][n] = gsl_cdf_beta_Qinv(beta_conf, s + 1/2, n - s + 1/2).
   This is the high bound of the Jeffrey's posterior for binomial
   confidence intervals. */
double beta_hi[NUM_BETA_PRECALC][NUM_BETA_PRECALC];


void init_beta(double beta_conf)
{
    int f, n;
    double ds, dn, beta_conf_inv = 1.0 - beta_conf;
    for (n = 0; n != NUM_BETA_PRECALC; ++n)
        for (f = 0; f <= n; ++f)
        {
            ds = (double)(n - f);
            dn = (double)n;
            beta_lo[f][n] = beta_Pinv(beta_conf_inv, ds + 0.5, dn - ds + 0.5);
            beta_hi[f][n] = beta_Qinv(beta_conf_inv, ds + 0.5, dn - ds + 0.5);
        }
}


/* safe function for obtaining a beta value */
static inline double jeffreys_beta_lo(int n, int s, double beta_conf)
{
    return n < NUM_BETA_PRECALC
        ? beta_lo[n - s][n]
        : beta_Pinv(1 - beta_conf, (double)s + 0.5, (double)(n - s) + 0.5);
}

static inline double jeffreys_beta_hi(int n, int s, double beta_conf)
{
    return n < NUM_BETA_PRECALC
        ? beta_hi[n - s][n]
        : beta_Qinv(1 - beta_conf, (double)s + 0.5, (double)(n - s) + 0.5);
}


/* mnemonics for labeling the low and high bounds for the beta
   estimate. */
enum bound_class {
    BOUND_CHANGED,
    BOUND_AMBIGUOUS,
    BOUND_UNCHANGED
};


/* Sample pairs of points from dist_pair up to max_points, classifying
   each pair as 'success' if distance is less than min_dist, 'failure'
   otherwise.  From the set of successes and failures, use the Beta
   distribution to estimate the true binomial probability.  Use pgen1
   and pgen2 to generate more points (and weights) as needed. */
enum fuzzy_state
binomial_quantile_est(unsigned max_points, float min_dist,
                      float post_conf, float beta_conf,
                      struct points_gen pgen1,
                      struct points_buf *points1,
                      struct points_gen pgen2,
                      struct points_buf *points2,
                      size_t batch_size)
{
    
    int n = 0, s = 0; /* # samples taken, # successes */

    enum bound_class lo_tag = BOUND_CHANGED, hi_tag = BOUND_UNCHANGED;

    /* min and max quantiles for posterior */
    float post_qmin = 1.0 - post_conf, post_qmax = post_conf;

    /* possibly weight-adjusted quantile values corresponding to
       beta_qmin and beta_qmax */
    float beta_qvmin = 0, beta_qvmax = 1;
    float dist_squared, min_dist_squared = gsl_pow_2(min_dist);

    assert(max_points % batch_size == 0);
    assert(points1->alloc >= max_points);
    assert(points2->alloc >= max_points);

    /* cur: next point to be used for distance calculation.
       end: next point to be drawn from distribution */
    POINT 
        *pcur1 = points1->buf + n,
        *pcur2 = points2->buf + n,
        *pend1 = points1->buf + points1->size,
        *pend2 = points2->buf + points2->size;
    /* invariant: pcur <= pend */
    
    unsigned p;
    while (n != max_points && (lo_tag != hi_tag || (beta_qvmax - beta_qvmin) > 0.05))
    {
        /* process another batch of samples, generating sample
           points as needed. */
        n += batch_size;
        while (points1->size < n) 
        {
            pgen1.gen_point(pgen1.gen_point_par, pend1);
            points1->size += batch_size;
            pend1 += batch_size;
        }
        while (points2->size < n) 
        {
            pgen2.gen_point(pgen2.gen_point_par, pend2);
            points2->size += batch_size;
            pend2 += batch_size;
        }
        
        /* measure distances, threshold, and classify successes */
        for (p = 0; p != batch_size; ++p)
        {
            dist_squared = 0;
            int d;
            for (d = 0; d != NUM_NUCS; ++d)
                dist_squared += gsl_pow_2((*pcur1)[d] - (*pcur2)[d]);

            s += (dist_squared < min_dist_squared ? 1 : 0);
            ++pcur1;
            ++pcur2;
        }

        /* regardless of state, we need to calculate beta_qvmin */
        beta_qvmin = jeffreys_beta_lo(n, s, beta_conf);
        
        if (post_qmax < beta_qvmin) lo_tag = BOUND_UNCHANGED;
        else if (post_qmin < beta_qvmin) lo_tag = BOUND_AMBIGUOUS;
        else lo_tag = BOUND_CHANGED;

        if (lo_tag != BOUND_UNCHANGED)
        {
            /* Now, calculate beta_qvmax only if necessary */
            beta_qvmax = jeffreys_beta_hi(n, s, beta_conf);
            assert(!isnan(beta_qvmax));

            if (post_qmax < beta_qvmax) hi_tag = BOUND_UNCHANGED;
            else if (post_qmin < beta_qvmax) hi_tag = BOUND_AMBIGUOUS;
            else hi_tag = BOUND_CHANGED;
        }
    }

    enum fuzzy_state state = AMBIGUOUS;

    if (lo_tag == hi_tag)
        switch(lo_tag)
        {
        case BOUND_CHANGED: state = CHANGED; break;
        case BOUND_AMBIGUOUS: state = AMBIGUOUS; break;
        case BOUND_UNCHANGED: state = UNCHANGED; break;
        }
    else if (lo_tag < hi_tag)
        state = (lo_tag == BOUND_CHANGED)
            ? AMBIGUOUS_OR_CHANGED
            : AMBIGUOUS_OR_UNCHANGED;
    else
    {
        fprintf(stderr, "%s:%i: low and hi bounds cross each other\n", __FILE__, __LINE__);
        exit(1);
    }
    return state;
}


#if 0
/* using existing points, calculate weights as necessary, apply
   threshold to classify as success or failure, and compute the
   weight-based modification to the Jeffrey's error estimates */
void weighted_binomial_est(struct points_gen pgen1,
                           struct weights_buf *weights1,
                           struct points_gen pgen2,
                           struct weights_buf *weights2)
{
    float w_succ = 0, w_fail = 0; /* sum of weights of successes and failures */
    float c_succ, c_fail; /* correction factors */

    /* cur: next weight to be used for distance calculation.
       end: next weight to be calculated from distribution ratios. */
    double
        *wcur1 = weights1->buf + n,
        *wcur2 = weights2->buf + n,
        *wend1 = weights1->buf + weights1->size,
        *wend2 = weights2->buf + weights2->size;

    double ww; /* product of weights */


    while (weights1->size < n)
    {
        pgen1.weight(points1->buf + weights1->size, pgen1.weight_par, wend1);
        weights1->size += batch_size;
        wend1 += batch_size;
    }
    while (weights2->size < n)
    {
        pgen2.weight(points2->buf + weights2->size, pgen2.weight_par, wend2);
        weights2->size += batch_size;
        wend2 += batch_size;
    }

    w = wcur1 - weights1->buf;
    while (wcur1 != weights1->buf + n)
    {
        ww = *wcur1++ * *wcur2++;
        if (success[w++]) w_succ += ww;
        else w_fail += ww;
    }

    c_succ = w_succ / (float)s;
    c_fail = w_fail / (float)(n - s);

    /* adjust beta_qvmin using the weights */
    beta_qvmin = (c_fail * beta_qvmin) 
        / (c_succ * (1.0 - beta_qvmin) + (c_fail * beta_qvmin));

    if (s != 0 && s != n)
        beta_qvmax = (c_fail * beta_qvmax)
            / (c_succ * (1.0 - beta_qvmax) + (c_fail * beta_qvmax));

}
#endif
