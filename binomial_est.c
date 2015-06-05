#include "binomial_est.h"

#include "geometry.h"

#include <gsl/gsl_cdf.h>
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_errno.h>

#include <stdlib.h>
#include <assert.h>
#include <stdio.h>
#include <pthread.h>

#include "cache.h"
#include "dirichlet_points_gen.h"

const char *fuzzy_state_strings[] = {
    "changed",
    "ambiguous_or_changed",
    "ambiguous",
    "ambiguous_or_unchanged",
    "unchanged"
};

/*
See the section 'Jeffreys Interval' on estimation of confidence
intervals on p given a sequence of bernoulli
outcomes. 

http://en.wikipedia.org/wiki/Binomial_proportion_confidence_interval


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

/* From gsl betainv.c. 
   Using binary seach, find x in [x0, x1] such that:
   abs(gsl_cdf_beta_P(x, a, b) - P) < Ptol
   or [x0, x1] shrinks to less than xtol.
 */
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
    x = bisect(x, P, a, b, 1e-10, 1e-10);
    return x;

}


double beta_Qinv(double Q, double a, double b)
{
    if (Q > 0.5) return beta_Pinv(1 - Q, a, b);
    else return 1 - beta_Pinv(Q, b, a);
}


/* Maximum N of precalculated Beta values. */


/* beta_lo[f][b] = gsl_cdf_beta_Pinv(beta_conf, s + 1/2, n - s + 1/2),
   where s: # of successes, f: # of failures, (n - s) n: # of trials,
   b = n / batch_size. This is the low bound of the Jeffrey's
   posterior for binomial confidence intervals. */

/* beta_hi[f][b] = gsl_cdf_beta_Qinv(beta_conf, s + 1/2,
   n - s + 1/2).  This is the high bound of the Jeffrey's posterior
   for binomial confidence intervals. */
struct {
    float *lo_buf, **lo, *hi_buf, **hi;
    unsigned batch_sz, max_n;
    double conf;
} beta_cache;


struct beta_input {
    unsigned start_batch, jump;
};

void *init_beta_func(void *args)
{
    struct beta_input *bi = args;
    double Q = 1.0 - beta_cache.conf;
    assert(Q < 0.5);
    unsigned f, n, b, n_batch = beta_cache.max_n / beta_cache.batch_sz;
    
    for (b = bi->start_batch; b < n_batch; b += bi->jump)
    {
        n = b * beta_cache.batch_sz;
        for (f = 0; f <= n; ++f)
            beta_cache.lo[f][b] = 
                beta_Pinv(Q, (double)(n - f) + 0.5, (double)f + 0.5);

        for (f = 0; f <= n; ++f)
            beta_cache.hi[f][b] = 1.0 - beta_cache.lo[n - f][b];
        
    }
    pthread_exit(NULL);
}

#define CHECK_THREAD(t, rc)                                             \
    if (rc)                                                             \
    {                                                                   \
        fprintf(stderr,                                                 \
                "Couldn't create thread %u at %s:%u with return code %i\n", \
                t, __FILE__, __LINE__, rc);                             \
        exit(1);                                                        \
    }                                                                   \
    
void binomial_est_init(double beta_conf, 
                       unsigned batch_size, 
                       unsigned num_beta_precalc,
                       size_t n_threads)
{
    pthread_t *threads = malloc(n_threads * sizeof(pthread_t));
    struct beta_input *inputs = malloc(n_threads * sizeof(struct beta_input));

    beta_cache.max_n = num_beta_precalc;
    beta_cache.batch_sz = batch_size;
    beta_cache.conf = beta_conf;

    unsigned n_batch = beta_cache.max_n / beta_cache.batch_sz;

    beta_cache.lo_buf = malloc(beta_cache.max_n * n_batch * sizeof(beta_cache.lo_buf[0]));
    beta_cache.lo = malloc(beta_cache.max_n * sizeof(beta_cache.lo_buf));

    beta_cache.hi_buf = malloc(beta_cache.max_n * n_batch * sizeof(beta_cache.hi_buf[0]));
    beta_cache.hi = malloc(beta_cache.max_n * sizeof(beta_cache.hi_buf));

    float **p, **pe, *b;
    pe = beta_cache.lo + beta_cache.max_n;
    for (p = beta_cache.lo, b = beta_cache.lo_buf; p != pe; ++p, b += n_batch)
        *p = b;

    pe = beta_cache.hi + beta_cache.max_n;
    for (p = beta_cache.hi, b = beta_cache.hi_buf; p != pe; ++p, b += n_batch)
        *p = b;


    unsigned t;
    int rc;
    for (t = 0; t != n_threads; ++t)
    {
        inputs[t] = (struct beta_input){ t, n_threads };

        rc = pthread_create(&threads[t], NULL, init_beta_func, &inputs[t]);
        CHECK_THREAD(t, rc);
    }
    for (t = 0; t != n_threads; ++t) {
        rc = pthread_join(threads[t], NULL);
        CHECK_THREAD(t, rc);
    }

    /* unsigned f, bi, n; */
    /* for (bi = 0; bi != n_batch; ++bi) */
    /* { */
    /*     n = bi * beta_cache.batch_sz; */
    /*     for (f = 0; f <= n; ++f) */
    /*         fprintf(stderr, "%u\t%u\t%7.5g\t%7.5g\n", f, n,  */
    /*                 beta_cache.lo[f][bi],  */
    /*                 beta_cache.hi[f][bi]); */
    /* } */

    free(threads);
    free(inputs);
}


void binomial_est_free()
{
    free(beta_cache.lo_buf);
    free(beta_cache.lo);
    free(beta_cache.hi_buf);
    free(beta_cache.hi);
}

/* safe function for obtaining a beta value */
static inline double jeffreys_beta_lo(int n, int s)
{
    assert(n % beta_cache.batch_sz == 0);
    return n < beta_cache.max_n
        ? beta_cache.lo[n - s][n / beta_cache.batch_sz]
        : beta_Pinv(1 - beta_cache.conf, (double)s + 0.5, (double)(n - s) + 0.5);
}

static inline double jeffreys_beta_hi(int n, int s)
{
    assert(n % beta_cache.batch_sz == 0);
    return n < beta_cache.max_n
        ? beta_cache.hi[n - s][n / beta_cache.batch_sz]
        : beta_Qinv(1 - beta_cache.conf, (double)s + 0.5, (double)(n - s) + 0.5);
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
   and pgen2 to generate more points (and weights) as needed. 
   
   Stopping criteria: 

   We want to take enough points to achieve a certain level of desired
   confidence (loose_spread). But, in the case our confidence interval
   straddles a cut-point, we demand extra confidence
   (tight_spread). It makes some sense for loose_spread to be the same
   as post_qmin on a vague intuitive level.  */
struct binomial_est_state
binomial_quantile_est(unsigned max_points, 
                      float min_dist,
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
    struct binomial_est_state est;
    est.beta_qval_lo = 0;
    est.beta_qval_hi = 1;
    double min_dist_squared = gsl_pow_2(min_dist);
    double *square_dist_buf = malloc(batch_size * sizeof(double));

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
    /*
    unsigned *alpha1_cts = ((struct points_gen_par *)pgen1.points_gen_par)->alpha_counts;
    unsigned *alpha2_cts = ((struct points_gen_par *)pgen2.points_gen_par)->alpha_counts;

    fprintf(stderr, "binomial_est_state: %u,%u,%u,%u\t%u,%u,%u,%u\n",
            alpha1_cts[0], alpha1_cts[1], alpha1_cts[2], alpha1_cts[3],
            alpha2_cts[0], alpha2_cts[1], alpha2_cts[2], alpha2_cts[3]);
    */

    /* See NOTE above on stopping criteria. */
    unsigned p;
    double
        beta_spread = est.beta_qval_hi - est.beta_qval_lo, 
        spread_loose_max = post_qmin,
        spread_tight_max = post_qmin * 0.2;
    while (n != max_points 
           && ((lo_tag != CHANGED && beta_spread > spread_loose_max)
               || (lo_tag == CHANGED && beta_spread > spread_tight_max)))
    {
        /* process another batch of samples, generating sample
           points as needed. */
        n += batch_size;
        while (points1->size < n) 
        {
            pgen1.gen_point(pgen1.points_gen_par, pend1);
            points1->size += batch_size;
            pend1 += batch_size;
        }
        while (points2->size < n) 
        {
            pgen2.gen_point(pgen2.points_gen_par, pend2);
            points2->size += batch_size;
            pend2 += batch_size;
        }
        
        /* measure distances, threshold, and classify successes.
           success means non-change */
        compute_square_dist((const double *)pcur1, 
                            (const double *)pcur2, 
                            batch_size, NUM_NUCS, square_dist_buf);

        for (p = 0; p != batch_size; ++p)
            s += (square_dist_buf[p] < min_dist_squared ? 1 : 0);

        /* regardless of est, we need to calculate est.beta_qval_lo */
        est.beta_qval_lo = jeffreys_beta_lo(n, s);
        
        if (post_qmax < est.beta_qval_lo) lo_tag = BOUND_UNCHANGED;
        else if (post_qmin < est.beta_qval_lo) lo_tag = BOUND_AMBIGUOUS;
        else lo_tag = BOUND_CHANGED;

        if (lo_tag != BOUND_UNCHANGED)
        {
            /* Now, calculate est.beta_qval_hi only if necessary */
            est.beta_qval_hi = jeffreys_beta_hi(n, s);
            assert(!isnan(est.beta_qval_hi));

            if (post_qmax < est.beta_qval_hi) hi_tag = BOUND_UNCHANGED;
            else if (post_qmin < est.beta_qval_hi) hi_tag = BOUND_AMBIGUOUS;
            else hi_tag = BOUND_CHANGED;
        }
        beta_spread = est.beta_qval_hi - est.beta_qval_lo;
    }

    est.state = AMBIGUOUS;

    if (lo_tag == hi_tag)
        switch(lo_tag)
        {
        case BOUND_CHANGED: est.state = CHANGED; break;
        case BOUND_AMBIGUOUS: est.state = AMBIGUOUS; break;
        case BOUND_UNCHANGED: est.state = UNCHANGED; break;
        }
    else if (lo_tag < hi_tag)
        est.state = (lo_tag == BOUND_CHANGED)
            ? AMBIGUOUS_OR_CHANGED
            : AMBIGUOUS_OR_UNCHANGED;
    else
    {
        fprintf(stderr, "%s:%i: low and hi bounds cross each other\n", __FILE__, __LINE__);
        exit(1);
    }
    free(square_dist_buf);
    return est;
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

    /* adjust est.beta_qval_lo using the weights */
    est.beta_qval_lo = (c_fail * est.beta_qval_lo) 
        / (c_succ * (1.0 - est.beta_qval_lo) + (c_fail * est.beta_qval_lo));

    if (s != 0 && s != n)
        est.beta_qval_hi = (c_fail * est.beta_qval_hi)
            / (c_succ * (1.0 - est.beta_qval_hi) + (c_fail * est.beta_qval_hi));

}
#endif
