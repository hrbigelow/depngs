#include "founder_dist.h"
/* A set of functions for calculating the distribution of the hidden
   founder base distribution from a set of basecalls and binned
   quality scores. */


/* number of quality score bins, configured at run-time. */
unsigned g_n_qual_bins = 8;
double g_qual_bin_error[MAX_QUAL_BINS];


/* compute the likelihood of hidden state composition with 'hs'
   secondary symbols, where the obsered configuration has op primary
   symbols and os secondary symbols, and the 'error' (miscall)
   probability is e. */
double
path_likelihood(unsigned op, unsigned os, unsigned hs, double e)
{
    /* strategy: sum in the direction of lengening path. factor out
       the initial binomial coefficients for the two components as
       well as the initial path contribution, yielding a first term of
       1. each successive term is updated from this first term using a
       ratio technique. the final sum is then adjusted by the initial
       factors. */
    /* # of secondary-to-primary symbols (hidden-to-observed) */
    int nsp = hs > op ? op : hs;

    /* # of secondary-to-secondary symbols */
    int nss = hs < os ? os - hs : os;

    double init_p_binom = m_choose_n(op, nsp);
    double init_s_binom = m_choose_n(os, nss);
    double path_factor = e / (1.0 - e); /* change in path prob in
                                           lengthing by one. */
    double path_factor_sq = path_factor * path_factor;
    double init_path_prob = 
        gsl_pow_int(e, nsp + (os - nss))
        * gsl_pow_int(1.0 - e, nss + (op - nsp));

    double marg = 0, term = 1.0;
    for (; nss >= 0 && nsp <= op; ++nsp, --nss) {
        /* replace with lookup table? */
        double p_binom_rat = (double)
        double s_binom_rat = (double)os / (double)(os - nss + 1.0);
        term *= p_binom_rat * s_binom_rat * path_factor_sq;
        marg += term;
    }
    return marg * init_path_prob * init_p_binom * init_s_binom;
}


/* find the mode of a component path function F(x;op;e) =
   m_choose_n(op,x) * e^(op-x)(1-e)^x using the ratio technique. Find
   the lowest value of x for which F(x+1)/F(x) < 1. */
unsigned
component_mode(unsigned op, double e)
{
    /* Derivation: F(x+1)/F(x) = {(x+1)/(op-x)} * {(1-e)/e}. Find the
       continuous value of x such that F(x+1)/F(x) = 1, then take the
       ceiling of x.  Let C = (1-e)/e. Solving for x, (op-x)/(x+1) =
       C. (op-x) = Cx+C. op-C = (C+1)x. (op-C)/(C+1) = x. */

    float C = (1.0-e)/e;
    return (unsigned)ceil((op-C)/(C+1));
}


struct cond_dist {
    unsigned maj_start, n;
    double dist[FLEX_ARRAY];
};

/* global parameter used for determining when to truncate a component
   cond distribution. if a given point carries less than
   g_frac_component_trunc fraction of the total distribution's mass in
   the majority set, it is truncated and not recorded. */
static double g_frac_component_trunc = 1e-4;

/* compute the majority set of the conditional distribution over the
   composition of H, Cond(H;op,os,e). */
struct cond_dist *
hidden_cond_dist(unsigned op, unsigned os, double e)
{
    /* strategy: partially fill a large buffer with the */
    double *dist_tmp = malloc(sizeof(double) * (os + op));

    /* find the component mode */
    unsigned hs = component_mode(op, e);
    unsigned end = hs + os;
    int s;
    unsigned p, s;
    double tot = 0;
    for (s = hs; s != end; ++s) {
        dist_tmp[s] = path_likelihood(op, os, s, e);
        tot += dist_tmp[s];
    }
    unsigned n = os;

    double comp, thresh = tot * g_frac_component_trunc; /* !!! parameter */

    /* search to the left and right */
    for (s = hs - 1; s >= 0, --s) {
        if ((comp = path_likelihood(op, os, s, e)) < thresh) break;
        dist_tmp[s] = comp;
        ++n;
    }
    unsigned start = s;
    for (s = end; s != op + os; ++s) {
        if ((comp = path_likelihood(op, os, s, e)) < thresh) break;
        dist_tmp[s] = comp;
        ++n;
    }
    struct cond_dist *dist = malloc(sizeof(struct cond_dist) + sizeof(double) * n);
    dist->maj_start = s;
    dist->n = n;
    memcpy(dist->dist, dist_tmp + s, sizeof(double) * n);
    free(dist_tmp);
    return dist;
}


KHASH_MAP_INIT_INT(hidden_dist_h, struct cond_dist *);
// typedef khash_t(hidden_dist_h) hidden_dist_t;

/* global, pre-populated hash of hidden-state conditional
   distributions. after the initial pre-computation phase, this hash
   will remain constant throughout the program run, in order to avoid
   mutexes. */
static hidden_dist_t g_hidden_dist_hash;


/* combine individual hidden-state conditional probability components
   into one distribution that describes the full conditional across
   the locus data. */
struct cond_dist *
combine_components(struct binned_qual_counts bqc)
{
    /* pre-allocate a cond_dist of the appropriate size */
    struct cond_dist *comps[MAX_QUAL_BINS];

    /* retrieve or calculate each individual component distribution */
    unsigned qb;
    double e; /* error probability */
    khiter_t itr;
    struct founder_pair_counts pc;
    struct cond_dist *comp;

    #define SCRATCH 100
    double dist_scratch[SCRATCH];

    /* copy the contents of the first */
    pc = bqc.bins[0];
    itr = kh_get(hidden_dist_h, g_hidden_dist_hash, pc);
    if (itr == kh_end(g_hidden_dist_hash)) {
        comp = hidden_cond_dist(pc.major, pc.minor, e);
        memcpy();
    }

    for (qb = 1; qb != g_n_qual_bins; ++qb) {
        pc = bqc.bins[qb];
        if (pc.major == 0 && pc.minor == 0) continue;
        e = g_qual_bin_error[qb]; /* lookup table */
        itr = kh_get(hidden_dist_h, g_hidden_dist_hash, pc);
        if (itr == kh_end(g_hidden_dist_hash)) {
            comp = hidden_cond_dist(pc.major, pc.minor, e);


        comps[qb] = 
    }

    struct cond_dist *dist = malloc(
    /* go through each individual component. */
}
