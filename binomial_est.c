#include "binomial_est.h"

#include <gsl/gsl_cdf.h>
#include <gsl/gsl_math.h>

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

/* mnemonics for labeling the low and high bounds for the beta
   estimate. */
enum bound_class {
    BOUND_CHANGED,
    BOUND_AMBIGUOUS,
    BOUND_UNCHANGED
};


/* Maximum N of precalculated Beta values. */
#define NUM_BETA_PRECALC 100

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
    double ds, dn;
    for (n = 0; n != NUM_BETA_PRECALC; ++n)
        for (f = 0; f != n; ++f)
        {
            ds = (double)(n - f);
            dn = (double)n;
            beta_lo[f][n] = gsl_cdf_beta_Pinv(beta_conf, ds + 0.5, dn - ds + 0.5);
            beta_hi[f][n] = gsl_cdf_beta_Qinv(beta_conf, ds + 0.5, dn - ds + 0.5);
        }
}


/* safe function for obtaining a beta value */
inline double jeffreys_beta_lo(int n, int s, double beta_conf)
{
    return n < NUM_BETA_PRECALC
        ? beta_lo[n - s][n]
        : gsl_cdf_beta_Pinv(beta_conf, (double)s + 0.5, (double)(n - s) + 0.5);
}

inline double jeffreys_beta_hi(int n, int s, double beta_conf)
{
    return n < NUM_BETA_PRECALC
        ? beta_hi[n - s][n]
        : gsl_cdf_beta_Qinv(beta_conf, (double)s + 0.5, (double)(n - s) + 0.5);
}

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
                      struct weights_buf *weights1,
                      struct points_gen pgen2,
                      struct points_buf *points2,
                      struct weights_buf *weights2,
                      size_t batch_size)
{
    
    int n = 0, s = 0; /* # samples taken, # successes */
    float w_succ = 0, w_fail = 0; /* sum of weights of successes and failures */
    float c_succ, c_fail; /* correction factors */

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
    assert(weights1->alloc >= max_points);
    assert(weights2->alloc >= max_points);

    POINT 
        *pcur1 = points1->buf + n,
        *pcur2 = points2->buf + n,
        *pend1 = points1->buf + points1->size,
        *pend2 = points2->buf + points2->size;

    double
        *wcur1 = weights1->buf + n,
        *wcur2 = weights2->buf + n,
        *wend1 = weights1->buf + weights1->size,
        *wend2 = weights2->buf + weights2->size;
    
    int p;

    unsigned char *success = (unsigned char *)malloc(batch_size);

    while (n != max_points && lo_tag != hi_tag)
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

            s += (success[p] = dist_squared < min_dist_squared ? 1 : 0);
            ++pcur1;
            ++pcur2;
        }

        /* regardless of state, we need to calculate beta_qvmin */
        beta_qvmin = jeffreys_beta_lo(n, s, beta_conf);

        if (s != 0 && s != n)
        {
            /* have at least one point on either side of cut, so
             weights can be informative. 
             !!! We might need to keep separate track of weights here...
            */
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

            for (p = 0; p != batch_size; ++p)
            {
                if (success[p]) w_succ += weights1[p] * weights2[p];
                else w_fail += weights1[p] * weights2[p];
            }
            c_succ = w_succ / (float)s;
            c_fail = w_fail / (float)(n - s);

            /* adjust beta_qvmin using the weights */
            beta_qvmin = (c_fail * beta_qvmin) 
                / (c_succ * (1.0 - beta_qvmin) + (c_fail * beta_qvmin));

        }
        
        if (post_qmax < beta_qvmin) lo_tag = BOUND_UNCHANGED;
        else if (post_qmin < beta_qvmin) lo_tag = BOUND_AMBIGUOUS;
        else lo_tag = BOUND_CHANGED;

        if (lo_tag != BOUND_UNCHANGED)
        {
            /* Now, calculate beta_qvmax only if necessary */
            beta_qvmax = jeffreys_beta_hi(n, s, beta_conf);
            
            if (s != 0 && s != n)
                beta_qvmax = (c_fail * beta_qvmax)
                    / (c_succ * (1.0 - beta_qvmax) + (c_fail * beta_qvmax));
            
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
        state = lo_tag == BOUND_CHANGED 
            ? AMBIGUOUS_OR_CHANGED
            : AMBIGUOUS_OR_UNCHANGED;
    else
    {
        fprintf(stderr, "%s:%i: low and hi bounds cross each other\n", __FILE__, __LINE__);
        exit(1);
    }
    free(success);
    return state;
}
