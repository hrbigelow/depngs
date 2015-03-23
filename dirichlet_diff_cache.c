/* provides memoization of classifying two dirichlets as same or
   different, based on their two largest components, and various
   confidence thresholds. */

#include "binomial_est.h"
#include "dist_worker.h"

#define MAX_COUNT1 1000
#define MAX_COUNT2 100

struct binomial_est_bounds {
    int16_t ambiguous[2];
    int16_t unchanged[2];
};

struct binomial_est_bounds bounds_cache[MAX_COUNT2][MAX_COUNT2][MAX_COUNT1];


struct pair_point_gen {
    struct points_gen pgen1, pgen2;
    struct points_buf pts1, pts2;
};

/* find top 2 components in el, store their indices in i1 and i2.
   Count the number of occurrences of min_val in z */
find_top2_aux(double *el, int *i1, int *i2, unsigned *z)
{
    double m1 = -1, m2 = -1;
    *i1 = -1, *i2 = -1;
    *z = 0;
    unsigned i;
    for (i = 0; i != NUM_NUCS; ++i)
    {
        if (m1 < el[i]) {
            *i2 = *i1, *i1 = i;
            m2 = m1, m1 = a[*i1];
        }
        if (el[i] == 0) ++*z;
    }
}

#define CACHE_GEN_POINTS_BATCH 20

/* Perform a binary search to locate the bounds.  Assume
   
 */
void initialize_est_bounds(unsigned b1, unsigned a2, unsigned b2,
                           double *prior_alpha, float min_dist,
                           float post_conf, float beta_conf,
                           unsigned max_points,
                           struct pair_point_gen *ppg,
                           struct binomial_est_bounds *beb)
{
    /* Do the binary search here. */
    enum fuzzy_state state;
    beb->ambiguous[0] = MAX_COUNT1 + 1;
    beb->ambiguous[1] = -1;
    beb->unchanged[0] = MAX_COUNT1 + 1;
    beb->unchanged[1] = -1;

    unsigned a1 = MAX_COUNT1 / 2;
    struct gen_dirichlet_points_par *gd1, *gd2;

    gd1 = ppg->pgen1.gen_point_par;
    memcpy(gd1->alpha, prior_alpha, sizeof(gd1->alpha));
    gd1->alpha[1] += a2;
    
    gd2 = ppg->pgen2.gen_point_par;
    memcpy(gd2->alpha, prior_alpha, sizeof(gd2->alpha));
    gd2->alpha[0] += b1;
    gd2->alpha[1] += b2;
    
    while (1)
    {
        gd1->alpha[0] = prior_alpha[0] + a1;
        state = 
            binomial_quantile_est(max_points, min_dist, post_conf, beta_conf,
                                  ppg->pgen1, ppg->pts1, ppg->pgen2, ppg->pts2,
                                  CACHE_GEN_POINTS_BATCH);
        if (state == AMBIGUOUS)
        {
            beb->ambiguous[0] = MIN(beb->ambiguous[0], a1);
            beb->ambiguous[1] = MAX(beb->ambiguous[1], a1);
        }
        else if (state == UNCHANGED)
        {
            beb->unchanged[0] = MIN(beb->unchanged
        }
}

/* test two dirichlets based on their counts. */
enum fuzzy_state cached_dirichlet_diff(unsigned *a_counts,
                                       unsigned *b_counts,
                                       double *prior_alpha,
                                       float min_dist,
                                       float post_conf,
                                       float beta_conf,
                                       unsigned max_points)
{
    int a1, a2, b1, b2;
    unsigned za, zb;
    find_top2_aux(a_counts, &a1, &a2, &za);
    find_top2_aux(b_counts, &b1, &b2, &zb);
    if (a1 == b1 && za >= 2 && zb >= 2
        && a_counts[a1] < MAX_COUNT1
        && a_counts[a2] < MAX_COUNT2
        && b_counts[b2] < MAX_COUNT2)
    {
        /* this qualifies for caching */
        struct binomial_est_bounds *beb = 
            &bounds_cache[a_counts[a2]][b_counts[b2]][b_counts[b1]];
        
        while (1)
        {
            while (beb->state == PENDING) ;
            if (beb->state == SET) return locate_cell(beb, a_counts[a1]);
            else
            {
                /* UNSET */
                pthread_mutex_lock(&set_flag_mtx);
                if (beb->state == PENDING)
                {
                    pthread_mutex_unlock(&set_flag_mtx);
                    continue;
                }
                else if (beb->state == SET)
                {
                    pthread_mutex_unlock(&set_flag_mtx);
                    return locate_cell(beb, a_counts[a1]);
                }
                else
                {
                    /* still UNSET */
                    beb->state = PENDING;
                    pthread_mutex_unlock(&set_flag_mtx);
                    initialize_est_bounds(a_counts[a1], a_counts[a2], b_counts[b2],
                                          prior_alpha, min_dist,
                                          post_conf, beta_conf, max_points, beb);
                    pthread_mutex_lock(&set_flag_mtx);
                    beb->state = SET;
                    pthread_mutex_unlock(&set_flag_mtx);
                    return locate_cell(beb, a_counts[a1]);
                }
            }
        }            
    }
    else
    {
        /* need to test the bounds organically, with no caching */
        
    }
}
