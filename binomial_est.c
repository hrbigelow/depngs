#include "binomial_est.h"

#include "geometry.h"

#include <gsl/gsl_cdf.h>
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_errno.h>

#include <assert.h>
#include <stdio.h>
#include <pthread.h>

#include "cache.h"

static struct binomial_est_params g_be_par;

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
static struct {
    float *lo_buf, **lo, *hi_buf, **hi;
    unsigned batch_sz, max_n;
    double conf;
} g_beta;


struct beta_input {
    unsigned start_batch, jump;
};

void *
init_beta_func(void *args)
{
    struct beta_input *bi = args;
    double Q = 1.0 - g_beta.conf;
    assert(Q < 0.5);
    unsigned f, n, b, n_batch = g_beta.max_n / g_beta.batch_sz;
    
    for (b = bi->start_batch; b < n_batch; b += bi->jump) {
        n = b * g_beta.batch_sz;
        for (f = 0; f <= n; ++f)
            g_beta.lo[f][b] = 
                beta_Pinv(Q, (double)(n - f) + 0.5, (double)f + 0.5);

        for (f = 0; f <= n; ++f)
            g_beta.hi[f][b] = 1.0 - g_beta.lo[n - f][b];
        
    }
    return NULL;
}

#define CHECK_THREAD(t, rc)                                             \
    if (rc)                                                             \
    {                                                                   \
        fprintf(stderr,                                                 \
                "Couldn't create thread %u at %s:%u with return code %i\n", \
                t, __FILE__, __LINE__, rc);                             \
        exit(1);                                                        \
    }                                                                   \
    
void
binomial_est_init(struct binomial_est_params be_par,
                  unsigned num_beta_precalc,
                  size_t n_threads)
{
    g_be_par = be_par;
    pthread_t *threads = malloc(n_threads * sizeof(pthread_t));
    struct beta_input *inputs = malloc(n_threads * sizeof(struct beta_input));

    g_beta.max_n = num_beta_precalc;
    g_beta.batch_sz = g_be_par.batch_size;
    g_beta.conf = g_be_par.beta_confidence;

    unsigned n_batch = g_beta.max_n / g_beta.batch_sz;

    g_beta.lo_buf = malloc(g_beta.max_n * n_batch * sizeof(g_beta.lo_buf[0]));
    g_beta.lo = malloc(g_beta.max_n * sizeof(g_beta.lo_buf));

    g_beta.hi_buf = malloc(g_beta.max_n * n_batch * sizeof(g_beta.hi_buf[0]));
    g_beta.hi = malloc(g_beta.max_n * sizeof(g_beta.hi_buf));

    float **p, **pe, *b;
    pe = g_beta.lo + g_beta.max_n;
    for (p = g_beta.lo, b = g_beta.lo_buf; p != pe; ++p, b += n_batch)
        *p = b;

    pe = g_beta.hi + g_beta.max_n;
    for (p = g_beta.hi, b = g_beta.hi_buf; p != pe; ++p, b += n_batch)
        *p = b;


    unsigned t;
    int rc;
    for (t = 0; t != n_threads; ++t) {
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
    /*     n = bi * g_beta.batch_sz; */
    /*     for (f = 0; f <= n; ++f) */
    /*         fprintf(stderr, "%u\t%u\t%7.5g\t%7.5g\n", f, n,  */
    /*                 g_beta.lo[f][bi],  */
    /*                 g_beta.hi[f][bi]); */
    /* } */

    free(threads);
    free(inputs);
}


void
binomial_est_free()
{
    free(g_beta.lo_buf);
    free(g_beta.lo);
    free(g_beta.hi_buf);
    free(g_beta.hi);
}

/* safe function for obtaining a beta value */
static inline double
jeffreys_beta_lo(int n, int s)
{
    assert(n % g_beta.batch_sz == 0);
    return n < g_beta.max_n
        ? g_beta.lo[n - s][n / g_beta.batch_sz]
        : beta_Pinv(1 - g_beta.conf, (double)s + 0.5, (double)(n - s) + 0.5);
}

static inline double
jeffreys_beta_hi(int n, int s)
{
    assert(n % g_beta.batch_sz == 0);
    return n < g_beta.max_n
        ? g_beta.hi[n - s][n / g_beta.batch_sz]
        : beta_Qinv(1 - g_beta.conf, (double)s + 0.5, (double)(n - s) + 0.5);
}


static inline double
beta_mean(int n, int s)
{
    return (s + 0.5) / (n + 0.5);
}


/* mnemonics for labeling the low and high bounds for the beta
   estimate. */


/* set est->lo_bound_tag, est->hi_bound_tag, and est->state based on
   the posterior confidence and beta_qval_{lo,hi} */
void
set_fuzzy_state(struct binomial_est_state *est, float post_conf)
{
    float post_qmax = post_conf, post_qmin = 1.0 - post_conf;
    est->lo_bound_tag =
        post_qmax < est->beta_lo
        ? BOUND_UNCHANGED
        : (post_qmin < est->beta_lo
           ? BOUND_AMBIGUOUS
           : BOUND_CHANGED);
           
    est->hi_bound_tag =
        post_qmax < est->beta_hi
        ? BOUND_UNCHANGED
        : (post_qmin < est->beta_hi
           ? BOUND_AMBIGUOUS
           : BOUND_CHANGED);
           
    if (est->lo_bound_tag == est->hi_bound_tag)
        switch(est->lo_bound_tag) {
        case BOUND_CHANGED: est->state = CHANGED; break;
        case BOUND_AMBIGUOUS: est->state = AMBIGUOUS; break;
        case BOUND_UNCHANGED: est->state = UNCHANGED; break;
        }
    else if (est->lo_bound_tag < est->hi_bound_tag)
        est->state = (est->lo_bound_tag == BOUND_CHANGED)
            ? AMBIGUOUS_OR_CHANGED
            : AMBIGUOUS_OR_UNCHANGED;
    else {
        fprintf(stderr, "%s:%i: low and hi bounds cross each other\n", __FILE__, __LINE__);
        exit(1);
    }
}


/* generate n_add more bernoulli trials.  buf1 and buf2 contain n_add
   points generated from the two dirichlet distributions of
   interest. */
void
add_binomial_trials(POINT *buf1, POINT *buf2, unsigned n_add,
                    struct binomial_est_state *est)
{
    double default_buf[32];
    double *square_dist = n_add > 32
        ? malloc(n_add * sizeof(double))
        : default_buf;

    double min_dist_squared = gsl_pow_2(g_be_par.min_dirichlet_dist);

    unsigned s = est->n_success, n = est->n_trials;
    n += n_add;

    /* measure distances, threshold, and classify successes.
       NOTE:  success means non-change */
    compute_square_dist((const double *)buf1,
                        (const double *)buf2, 
                        n_add, NUM_NUCS, square_dist);

    /* tallying number of 'successes' (pairs of points that are
       below the distance threshold, i.e. non-changed) */
    unsigned p;
    for (p = 0; p != n_add; ++p)
        s += (square_dist[p] < min_dist_squared ? 1 : 0);
    
    /* regardless of est, we need to calculate est->beta_lo */
    est->beta_lo = jeffreys_beta_lo(n, s);
    est->beta_hi = jeffreys_beta_hi(n, s);
    est->beta_mean = beta_mean(n, s);
    est->n_trials = n;
    est->n_success = s;
    set_fuzzy_state(est, g_be_par.post_confidence);

    if (square_dist != default_buf)
        free(square_dist);
}


/* Sample pairs of points from dp1 and dp2 up to the maximum number of
   points, generating new batches of points as needed.  measure
   euclidean distance of each pair of points, and threshold to get
   'success' if distance is less than min_dist, 'failure' otherwise.

   From the set of successes and failures, use the Beta distribution
   to estimate the true binomial probability.
   
   Stopping criteria: 

   We want to take enough points to achieve a certain level of desired
   confidence (loose_spread). But, in the case our confidence interval
   straddles a cut-point, we demand extra confidence
   (tight_spread). dp1 and dp2 are assumed initialized.  this call
   will add additional points to dp1 and dp2 until  */
struct binomial_est_state
binomial_quantile_est(struct dir_points *dp1,
                      struct dir_points *dp2,
                      unsigned batch_sz)
{
    struct binomial_est_state est = {
        .beta_lo = 0,
        .beta_hi = 1,
        .n_trials = 0,
        .n_success = 0
    };
    set_fuzzy_state(&est, g_be_par.post_confidence);
    
    double
        beta_spread = est.beta_hi - est.beta_lo,
        spread_loose_max = 1.0 - g_be_par.post_confidence,
        spread_tight_max = spread_loose_max * 0.2;

    POINT *buf1 = dp1->data, *buf2 = dp2->data;
    while (est.n_trials != g_be_par.max_sample_points 
           && ((est.lo_bound_tag != BOUND_CHANGED && beta_spread > spread_loose_max)
               || (est.lo_bound_tag == BOUND_CHANGED && beta_spread > spread_tight_max)))
    {
        if (dp1->n_points < est.n_trials + batch_sz) {
            gen_dir_points(dp1->perm_alpha, buf1, batch_sz);
            dp1->n_points += batch_sz;
        }
        if (dp2->n_points < est.n_trials + batch_sz) {
            gen_dir_points(dp2->perm_alpha, buf2, batch_sz);
            dp2->n_points += batch_sz;
        }
        add_binomial_trials(buf1, buf2, batch_sz, &est);
        beta_spread = est.beta_hi - est.beta_lo;
        buf1 += batch_sz;
        buf2 += batch_sz;
    }
    return est;
}
 


/* Interpolate the interval in the beb row */
enum fuzzy_state
locate_cell(struct binomial_est_bounds *beb, unsigned a1)
{
    if (beb->unchanged[0] <= a1 && a1 < beb->unchanged[1])
        return UNCHANGED;
    else if (beb->ambiguous[0] <= a1 && a1 < beb->ambiguous[1])
        return AMBIGUOUS;
    else return CHANGED;
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
        weights1->size += g_be_par.batch_size;
        wend1 += g_be_par.batch_size;
    }
    while (weights2->size < n)
    {
        pgen2.weight(points2->buf + weights2->size, pgen2.weight_par, wend2);
        weights2->size += g_be_par.batch_size;
        wend2 += g_be_par.batch_size;
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

    /* adjust est.beta_lo using the weights */
    est.beta_lo = (c_fail * est.beta_lo) 
        / (c_succ * (1.0 - est.beta_lo) + (c_fail * est.beta_lo));

    if (s != 0 && s != n)
        est.beta_hi = (c_fail * est.beta_hi)
            / (c_succ * (1.0 - est.beta_hi) + (c_fail * est.beta_hi));

}


struct binomial_est_state
binomial_quantile_est(struct dir_points *dp1,
                      struct dir_points *dp2,
                      unsigned batch_sz)
{
    int n = 0, s = 0; /* # samples taken, # successes */

    enum bound_class lo_tag = BOUND_CHANGED, hi_tag = BOUND_UNCHANGED;
    
    /* min and max quantiles for posterior */
    float
        post_qmin = 1.0 - g_be_par.post_confidence, 
        post_qmax = g_be_par.post_confidence;

    /* possibly weight-adjusted quantile values corresponding to
       beta_qmin and beta_qmax */
    struct binomial_est_state est;
    est.beta_lo = 0;
    est.beta_hi = 1;
    double min_dist_squared = gsl_pow_2(g_be_par.min_dirichlet_dist);
    double *square_dist_buf = malloc(batch_sz * sizeof(double));

    assert(g_be_par.max_sample_points % batch_sz == 0);

    /* cur: next point to be used for distance calculation.
       end: next point to be drawn from distribution */
    POINT 
        *pcur1 = dp1->data, *pend1 = pcur1 + dp1->n_points,
        *pcur2 = dp2->data, *pend2 = pcur2 + dp2->n_points;

    assert(pcur1 != NULL);
    assert(pcur2 != NULL);

    /* See NOTE above on stopping criteria. */
    unsigned p;
    double
        beta_spread = est.beta_hi - est.beta_lo, 
        spread_loose_max = post_qmin,
        spread_tight_max = post_qmin * 0.2;
    while (n != g_be_par.max_sample_points 
           && ((lo_tag != BOUND_CHANGED && beta_spread > spread_loose_max)
               || (lo_tag == BOUND_CHANGED && beta_spread > spread_tight_max)))
    {
        /* process another batch of samples, generating sample
           points as needed. */
        n += batch_sz;
        while (dp1->n_points < n) {
            gen_dir_points(dp1->perm_alpha, pend1, batch_sz);
            dp1->n_points += batch_sz;
            pend1 += batch_sz;
        }
        while (dp2->n_points < n) {
            gen_dir_points(dp2->perm_alpha, pend2, batch_sz);
            dp2->n_points += batch_sz;
            pend2 += batch_sz;
        }
        
        /* measure distances, threshold, and classify successes.
           NOTE:  success means non-change */
        compute_square_dist((const double *)pcur1, 
                            (const double *)pcur2, 
                            batch_sz, NUM_NUCS, square_dist_buf);
        pcur1 += batch_sz;
        pcur2 += batch_sz;

        /* tallying number of 'successes' (pairs of points that are
           below the distance threshold, i.e. non-changed) */
        for (p = 0; p != batch_sz; ++p)
            s += (square_dist_buf[p] < min_dist_squared ? 1 : 0);

        /* regardless of est, we need to calculate est.beta_lo */
        est.beta_lo = jeffreys_beta_lo(n, s);
        
        if (post_qmax < est.beta_lo) lo_tag = BOUND_UNCHANGED;
        else if (post_qmin < est.beta_lo) lo_tag = BOUND_AMBIGUOUS;
        else lo_tag = BOUND_CHANGED;

        if (lo_tag != BOUND_UNCHANGED) {
            /* Now, calculate est.beta_hi only if necessary */
            est.beta_hi = jeffreys_beta_hi(n, s);
            assert(!isnan(est.beta_hi));

            if (post_qmax < est.beta_hi) hi_tag = BOUND_UNCHANGED;
            else if (post_qmin < est.beta_hi) hi_tag = BOUND_AMBIGUOUS;
            else hi_tag = BOUND_CHANGED;
        }
        beta_spread = est.beta_hi - est.beta_lo;
    }

    est.state = AMBIGUOUS;

    if (lo_tag == hi_tag)
        switch(lo_tag) {
        case BOUND_CHANGED: est.state = CHANGED; break;
        case BOUND_AMBIGUOUS: est.state = AMBIGUOUS; break;
        case BOUND_UNCHANGED: est.state = UNCHANGED; break;
        }
    else if (lo_tag < hi_tag)
        est.state = (lo_tag == BOUND_CHANGED)
            ? AMBIGUOUS_OR_CHANGED
            : AMBIGUOUS_OR_UNCHANGED;
    else {
        fprintf(stderr, "%s:%i: low and hi bounds cross each other\n", __FILE__, __LINE__);
        exit(1);
    }
    free(square_dist_buf);
    est.n_trials = n;
    est.n_success = s;
    return est;
}


#endif
