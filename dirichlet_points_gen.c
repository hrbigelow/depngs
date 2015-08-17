#include "dirichlet_points_gen.h"

#include <gsl/gsl_randist.h>
#include <gsl/gsl_sf_exp.h>

#include <string.h>
#include <math.h>
#include "yepMath.h"
#include "yepCore.h"
#include <assert.h>
#include <float.h>

#define MIN(a, b) ((a) < (b) ? (a) : (b))
#define MAX(a, b) ((a) < (b) ? (b) : (a))

static double alpha_prior;

static struct phred {
    float prob_right; /* probability that the true base matches the basecall */
    float prob_wrong; /* probability that the true base is one of the
                         three non-basecall bases. */
} error_probability[255];

/* used for */
static double g_log_dbl_max, g_log_dbl_min, g_log_dbl_range;

void
dirichlet_points_gen_init(double _alpha_prior)
{
    alpha_prior = _alpha_prior;

    unsigned q;
    double ep;
    for (q = 0; q != 255; ++q) {
        ep = pow(10.0, -(double)q / 10.0);
        error_probability[q] = (struct phred){ 1.0 - ep, ep / 3.0 };
    }
    g_log_dbl_max = log(DBL_MAX);
    g_log_dbl_min = log(DBL_MIN);
    g_log_dbl_range = g_log_dbl_max - g_log_dbl_min;
}

double get_alpha_prior()
{
    return alpha_prior;
}

/* compute GEN_POINTS_BATCH of unnormalized dirichlet values using the
   same alpha.  assume NUM_NUCS dimensions */
void ran_dirichlet_lnpdf_unnormalized(double *alpha, double *points, double *lndir)
{
    POINT
        lnpoints[GEN_POINTS_BATCH], 
        *lnp, 
        *lnpe = lnpoints + GEN_POINTS_BATCH,
        alpha_minus_one;
    
    unsigned i;
    for (i = 0; i != NUM_NUCS; ++i)
        alpha_minus_one[i] = alpha[i] - 1.0;

    /* treating lnpoints as individual double components */
    (void)yepMath_Log_V64f_V64f(points, (double *)lnpoints, 
                                sizeof(lnpoints) / sizeof(double));

    lnp = lnpoints;
    for (lnp = lnpoints; lnp != lnpe; ++lnp, ++lndir) {
        *lndir = alpha_minus_one[0] * (*lnp)[0]
            + alpha_minus_one[1] * (*lnp)[1]
            + alpha_minus_one[2] * (*lnp)[2]
            + alpha_minus_one[3] * (*lnp)[3];
    }
}


/* the following four functions can be used with binomial_est's struct
   points_gen. */
void gen_dirichlet_points_wrapper(const void *par, POINT *points)
{
    int i;
    struct points_gen_par *gd = par;
    double alpha[] = { 
        gd->alpha_counts[0] + alpha_prior,
        gd->alpha_counts[1] + alpha_prior,
        gd->alpha_counts[2] + alpha_prior,
        gd->alpha_counts[3] + alpha_prior 
    };

    for (i = 0; i != GEN_POINTS_BATCH; ++i) {
        gsl_ran_dirichlet(gd->randgen, NUM_NUCS, alpha, *points);
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
    char refbase = *(char *)par;
    static char nucs[] = "ACGT";
    int ref_ind = index(nucs, refbase) - nucs;
    
    int i;
    for (i = 0; i != GEN_POINTS_BATCH; ++i, ++points)
        memcpy(*points, ref_points[ref_ind], sizeof(ref_points[0]));
}


/* calculates log(likelihood) - log(dirichlet).  In this, the
   dirichlet prior is a common factor and so cancels.  */
void
calc_post_to_dir_logratio(POINT *points, const void *par, double *weights)
{
    int i;
    const struct points_gen_par *pd = par;
    const struct bqs_count *trm = pd->observed, *trm_end = trm + pd->n_observed;

    double 
        dotp[GEN_POINTS_BATCH], ldotp[GEN_POINTS_BATCH],
        llh[GEN_POINTS_BATCH], ldir[GEN_POINTS_BATCH];

    memset(llh, 0, sizeof(llh));

    POINT *p, *pe = points + GEN_POINTS_BATCH;

    /* llh will not include the prior Dirichlet.  Only the effect of
       the data itself. To correct for this, we must use a
       'residual_alpha' for the equivalent, perfect-quality Dirichlet
       proposal.*/
    struct phred ph;
    unsigned base_code;
    while (trm != trm_end) {
        if (trm->qual < pd->min_base_quality) {
            ++trm;
            continue;
        }
        ph = error_probability[trm->qual];
        double prob[] = { ph.prob_wrong, ph.prob_wrong, ph.prob_wrong, ph.prob_wrong };
        base_code = seq_nt16_int[(int)trm->base];
        assert(base_code < 4);
        prob[base_code] = ph.prob_right;
        for (p = points, i = 0; p != pe; ++p, ++i)
            dotp[i] = 
                (*p)[0] * prob[0] + (*p)[1] * prob[1]
                + (*p)[2] * prob[2] + (*p)[3] * prob[3];

        (void)yepMath_Log_V64f_V64f(dotp, ldotp, GEN_POINTS_BATCH);
        if (trm->ct > 1) 
            (void)yepCore_Multiply_IV64fS64f_IV64f(ldotp, (double)trm->ct, GEN_POINTS_BATCH);
        (void)yepCore_Add_V64fV64f_V64f(llh, ldotp, llh, GEN_POINTS_BATCH);
        ++trm;
    }

    /* This is the residual Dirichlet correction. */
    POINT residual_alpha;
    for (i = 0; i != NUM_NUCS; ++i)
        residual_alpha[i] = (double)pd->alpha_counts[i] - alpha_prior + 1.0;
    ran_dirichlet_lnpdf_unnormalized(residual_alpha, (double *)points, ldir);

    /* for (i = 0; i != GEN_POINTS_BATCH; ++i) */
    /*     fprintf(stderr, "%7.5g\t%7.5g\t%7.5g\t%7.5g\t%7.5g\t%7.5g\n", */
    /*             points[i][0], points[i][1], points[i][2], points[i][3], llh[i], ldir[i]); */
    enum YepStatus stat =
        yepCore_Subtract_V64fV64f_V64f(llh, ldir, weights, GEN_POINTS_BATCH);
    assert(stat == YepStatusOk);
}


/* exponentiate vals, scaling to avoid underflow or
   overflow. ultimately we want  */
void
batch_scaled_exponentiate(double *val, unsigned n_val)
{
    double min_val = DBL_MAX, max_val = -DBL_MAX;

    unsigned i;
    for (i = 0; i != n_val; ++i) {
        min_val = MIN(min_val, val[i]);
        max_val = MAX(max_val, val[i]);
    }
    double adj;
    /* double val_spread = max_val - min_val; */
    adj = max_val;
    /* if (val_spread < g_log_dbl_range) */
    /*     adj = (max_val + min_val) * 0.5; */
    /* else */
    /*     adj = max_val - (g_log_dbl_max - 10); */

    enum YepStatus stat;
    stat = yepCore_Subtract_V64fS64f_V64f(val, adj, val, n_val);
    assert(stat == YepStatusOk);

    stat = yepMath_Exp_V64f_V64f(val, val, n_val);
    assert(stat == YepStatusOk);

    /* double min_val2 = DBL_MAX; */
    /* double max_val2 = -DBL_MAX; */
    /* for (i = 0; i != n_val; ++i) { */
    /*     min_val2 = MIN(min_val2, val[i]); */
    /*     max_val2 = MAX(max_val2, val[i]); */
    /* } */
    /* assert(max_val2 > 0); */
}


void calc_dummy_logratio(POINT *point, const void *par, double *weights)
{
    int i;
    for (i = 0; i != GEN_POINTS_BATCH; ++i)
        weights[i] = 0;
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
