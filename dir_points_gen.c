#include "dir_points_gen.h"

#include <gsl/gsl_randist.h>
#include <gsl/gsl_sf_exp.h>

#include <string.h>
#include <math.h>
#include <assert.h>
#include <float.h>

#include "dir_cache.h"
#include "geometry.h"

#include "yepLibrary.h"
#include "yepMath.h"
#include "yepCore.h"

#define MIN(a, b) ((a) < (b) ? (a) : (b))
#define MAX(a, b) ((a) < (b) ? (b) : (a))

static struct phred {
    float prob_right; /* probability that the true base matches the basecall */
    float prob_wrong; /* probability that the true base is one of the
                         three non-basecall bases. */
} error_probability[255];


/* used for */
static double g_log_dbl_max, g_log_dbl_min, g_log_dbl_range;

static struct dirichlet_points_gen_params g_pg_par;

static __thread gsl_rng *tls_rng;

void
dirichlet_points_gen_init(struct dirichlet_points_gen_params pg_par)
{
    g_pg_par = pg_par;
    enum YepStatus status = yepLibrary_Init();
    assert(status == YepStatusOk);

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


void
dirichlet_points_gen_free()
{
    enum YepStatus status = yepLibrary_Release();
    assert(status == YepStatusOk);
}


void
dir_points_thread_init()
{
    tls_rng = gsl_rng_alloc(gsl_rng_taus);
}


void
dir_points_thread_free()
{
    gsl_rng_free(tls_rng);
}


void
dir_points_update_alpha(const unsigned *alpha,
                        const unsigned *perm,
                        struct dir_points *dp)
{
    /* update the alphas and record whether there was a change */
    static unsigned perm_default[] = { 0, 1, 2, 3 };
    if (! perm) perm = perm_default;
    unsigned i, change = 0;
    for (i = 0; i != NUM_NUCS; ++i) {
        if (dp->perm_alpha[i] != alpha[perm[i]])
            change = 1;
        dp->perm_alpha[i] = alpha[perm[i]];
        dp->perm[i] = perm[i];
    }
    if (change) {
        dp->n_weights = 0;
        POINT *p = dir_cache_try_get_points(dp->perm_alpha);
        if (p) {
            dp->data = p;
            dp->n_points = g_pg_par.max_sample_points;
        } else {
            dp->data = dp->points_buf;
            dp->n_points = 0;
        }
    }
}


/* fully populate the points buffer if not full. call this after
   calling dir_points_update_alpha. */
void
dir_points_fill(struct dir_points *dp)
{
    unsigned msp = g_pg_par.max_sample_points;
    if (dp->n_points != msp) {
        gen_dir_points(dp->perm_alpha, dp->points_buf, msp);
        dp->data = dp->points_buf;
        dp->n_points = msp;
        dp->n_weights = 0;
    }
}


void
dir_weights_update_terms(struct bqs_count *bqs_ct, unsigned n_bqs_ct,
                         struct dir_points *dp)
{
    dp->bqs_ct = bqs_ct;
    dp->n_bqs_ct = n_bqs_ct;
}


/* must call de_permute_points before calling this. */
void
dir_weights_fill(struct dir_points *dp)
{
    assert(dp->perm[0] == 0 && dp->perm[1] == 1
           && dp->perm[2] == 2 && dp->perm[3] == 3);

    unsigned msp = g_pg_par.max_sample_points;
    assert(dp->n_points == msp);
    while (dp->n_weights != msp)
        calc_post_to_dir_logratio(dp);
    
    batch_scaled_exponentiate(dp->weights, dp->n_weights);
}


void
alloc_locus_data(struct locus_data *ld)
{
    ld->dist.points_buf = malloc(sizeof(ld->dist.points_buf[0]) * g_pg_par.max_sample_points);
    ld->dist.n_points = 0;
    ld->dist.data = NULL;
    ld->dist.weights = malloc(sizeof(ld->dist.weights[0]) * g_pg_par.max_sample_points);
    ld->dist.n_weights = 0;
    ld->bqs_ct = NULL;
    ld->indel_ct = NULL;
    init_pileup_data(&ld->sample_data);
}


void
free_locus_data(struct locus_data *ld)
{
    free(ld->dist.points_buf);
    free(ld->dist.weights);
    if (ld->bqs_ct != NULL) free(ld->bqs_ct);
    if (ld->indel_ct != NULL) free(ld->indel_ct);
    free_pileup_data(&ld->sample_data);
}


/* call when we advance to a new locus */
void
reset_locus_data(struct locus_data *ld)
{
    ld->init.base_ct = 0;
    ld->init.bqs_ct = 0;
    ld->init.indel_ct = 0;
    ld->init.sample_data = 0;
    ld->dist.n_points = 0;
    ld->dist.n_weights = 0;
    ld->confirmed_changed = 0;
    if (ld->bqs_ct != NULL) {
        free(ld->bqs_ct);
        ld->bqs_ct = NULL;
    }
    if (ld->indel_ct != NULL) {
        free(ld->indel_ct);
        ld->indel_ct = NULL;
    }
}


double get_alpha_prior()
{
    return g_pg_par.alpha_prior;
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


/* generate a complete set of dirichlet points */
void
gen_dir_points(unsigned *cts, POINT *points, unsigned n_points)
{
    double alpha[] = {
        cts[0] + g_pg_par.alpha_prior,
        cts[1] + g_pg_par.alpha_prior,
        cts[2] + g_pg_par.alpha_prior,
        cts[3] + g_pg_par.alpha_prior
    };
    while (n_points-- != 0) {
        gsl_ran_dirichlet(tls_rng, 4, alpha, (double *)points);
        ++points;
    }
}


/* re-order point components to the default permutation { 0, 1, 2, 3
   }.  handle the case where the dp is using the cache. */
void
de_permute_points(struct dir_points *dp)
{
    if (dp->perm[0] == 0 && dp->perm[1] == 1
        && dp->perm[2] == 2 && dp->perm[3] == 3)
        return;
    
    if (dp->data != dp->points_buf) {
        memcpy(dp->points_buf, dp->data, dp->n_points * sizeof(POINT));
        dp->data = dp->points_buf;
    }
    
    unsigned i, pinv[4];
    for (i = 0; i != 4; ++i)
        pinv[dp->perm[i]] = i;

    POINT *p, tmp;
    for (p = dp->data; p != dp->data + dp->n_points; ++p) {
        for (i = 0; i != 4; ++i) tmp[i] = (*p)[pinv[i]];
        for (i = 0; i != 4; ++i) (*p)[i] = tmp[i];
    }

    static unsigned perm_default[] = { 0, 1, 2, 3 };
    memcpy(dp->perm, perm_default, sizeof(perm_default));
}



/* calculates log(likelihood) - log(dirichlet).  In this, the
   dirichlet prior is a common factor and so cancels.  Since dp
   generally stores its points with coordinates in permuted form, we
   require that the points are de-permuted before this calculation. */
void
calc_post_to_dir_logratio(struct dir_points *dp)
{ 
    assert(dp->perm[0] == 0 && dp->perm[1] == 1
           && dp->perm[2] == 2 && dp->perm[3] == 3);

    int i;
    const struct bqs_count
        *term = dp->bqs_ct,
        *term_end = term + dp->n_bqs_ct;
    
    double 
        dotp[GEN_POINTS_BATCH], ldotp[GEN_POINTS_BATCH],
        llh[GEN_POINTS_BATCH], ldir[GEN_POINTS_BATCH];

    memset(llh, 0, sizeof(llh));
    POINT
        *p,
        *points = dp->data + dp->n_weights,
        *pe = points + GEN_POINTS_BATCH;

    /* start populating weights at the end */
    double *w = dp->weights + dp->n_weights;
    
    /* llh will not include the prior Dirichlet.  Only the effect of
       the data itself. To correct for this, we must use a
       'residual_alpha' for the equivalent, perfect-quality Dirichlet
       proposal.*/
    struct phred ph;
    unsigned base_code;
    while (term != term_end) {
        if (term->qual < g_pg_par.min_base_quality) {
            ++term;
            continue;
        }
        ph = error_probability[term->qual];
        double prob[] = { ph.prob_wrong, ph.prob_wrong, ph.prob_wrong, ph.prob_wrong };
        base_code = seq_nt16_int[(int)term->base];
        assert(base_code < 4);
        prob[base_code] = ph.prob_right;
        for (p = points, i = 0; p != pe; ++p, ++i)
            dotp[i] = 
                (*p)[0] * prob[0] + (*p)[1] * prob[1]
                + (*p)[2] * prob[2] + (*p)[3] * prob[3];

        (void)yepMath_Log_V64f_V64f(dotp, ldotp, GEN_POINTS_BATCH);
        if (term->ct > 1) 
            (void)yepCore_Multiply_IV64fS64f_IV64f(ldotp, (double)term->ct, GEN_POINTS_BATCH);
        (void)yepCore_Add_V64fV64f_V64f(llh, ldotp, llh, GEN_POINTS_BATCH);
        ++term;
    }

    /* This is the residual Dirichlet correction. */
    POINT residual_alpha;
    for (i = 0; i != NUM_NUCS; ++i)
        residual_alpha[i] = (double)dp->perm_alpha[i] - g_pg_par.alpha_prior + 1.0;
    ran_dirichlet_lnpdf_unnormalized(residual_alpha, (double *)points, ldir);

    /* for (i = 0; i != GEN_POINTS_BATCH; ++i) */
    /*     fprintf(stderr, "%7.5g\t%7.5g\t%7.5g\t%7.5g\t%7.5g\t%7.5g\n", */
    /*             points[i][0], points[i][1], points[i][2], points[i][3], llh[i], ldir[i]); */
    enum YepStatus stat =
        yepCore_Subtract_V64fV64f_V64f(llh, ldir, w, GEN_POINTS_BATCH);
    assert(stat == YepStatusOk);

    /* update the claim of number of valid weights */
    dp->n_weights += GEN_POINTS_BATCH;
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


/* populates square_dist_buf with squares of euclidean distances
   between points1 and points2 in the barycentric space (R4,
   normalized positive components).  populates weights_buf with
   product of weights1 and weights2. */
void
compute_wsq_dist(const double *points1, const double *weights1,
                 const double *points2, const double *weights2,
                 size_t n_points,
                 double *square_dist_buf, double *weights_buf)
{
    compute_square_dist(points1, points2, n_points, 4, square_dist_buf);
    (void)yepCore_Multiply_V64fV64f_V64f(weights1, weights2, weights_buf, n_points);
}
