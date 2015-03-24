/* provides memoization of classifying two dirichlets as same or
   different, based on their two largest components, and various
   confidence thresholds. */

#include "binomial_est.h"
#include "dist_worker.h"
#include "virtual_bound.h"


#define MAX_COUNT1 1000
#define MAX_COUNT2 100

/* bounds_cache[a2][b2][b1] = a description of the matrix of distance
   categories.  */
struct binomial_est_bounds bounds_cache[MAX_COUNT2][MAX_COUNT2][MAX_COUNT1];

/* Describes estimated distance between two Dirichlets with alphas equal to:
   { A1 + p, A2 + p, p, p }
   { B1 + p, B2 + p, p, p }

   Where A1, A2, B1, B2 are integral values.  [unchanged[0],
   unchanged[1]) represents the range of values of A1 where it is
   deemed UNCHANGED.  [ambiguous[0], ambiguous[1]) are the values of
   A1 for which it is deemed AMBIGUOUS (or UNCHANGED, where this
   interval overlaps the unchanged interval. */
struct binomial_est_bounds {
    int16_t ambiguous[2];
    int16_t unchanged[2];
};


/* Interpolate the interval in the beb row */
enum fuzzy_state locate_cell(struct binomial_est_bounds *beb, unsigned a1)
{
    if (beb->unchanged[0] <= a1 && a1 < beb->unchanged[1])
        return UNCHANGED;
    else if (beb->ambiguous[0] <= a1 && a1 < beb->ambiguous[1])
        return AMBIGUOUS;
    else return CHANGED;
}



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


int query_is_less(unsigned pos, void *par)
{
    struct binomial_est_params *b = par;
    ((struct gen_dirichlet_points_par *)b->pgen1->gen_points_par)->alpha[0] =
        (double)pos + b->prior_alpha[0];
    enum fuzzy_state elem_state = binomial_quantile_est(b->max_points, b->min_dist,
                                                        b->post_conf, b->beta_conf,
                                                        b->pgen1, b->points1,
                                                        b->pgen2, b->points2,
                                                        b->batch_size);
    return b->query_state < elem_state ? 1 : 0;
}


int elem_is_less(unsigned pos, void *par)
{
    struct binomial_est_params *b = par;
    ((struct gen_dirichlet_points_par *)b->pgen1->gen_points_par)->alpha[0] =
        (double)pos + b->prior_alpha[0];
    enum fuzzy_state elem_state = binomial_quantile_est(b->max_points, b->min_dist,
                                                        b->post_conf, b->beta_conf,
                                                        b->pgen1, b->points1,
                                                        b->pgen2, b->points2,
                                                        b->batch_size);
    return elem_state < b->query_state ? 1 : 0;
}




#define CACHE_GEN_POINTS_BATCH 20

void initialize_est_bounds(unsigned a2, unsigned b1, unsigned b2,
                           struct binomial_est_params *bpar,
                           struct binomial_est_bounds *beb)
{
    beb->ambiguous[0] = beb->ambiguous[1] = 0;
    beb->unchanged[0] = beb->unchanged[1] = 0;

    struct gen_dirichlet_points_par *gd1, *gd2;
    unsigned a1 = (unsigned)round(a2 + a2 * (b1 + b2) / (double)b2);

    gd1 = bpar->pgen1.gen_point_par;
    memcpy(gd1->alpha, bpar->prior_alpha, sizeof(gd1->alpha));
    gd1->alpha[0] += nearest_a1;
    gd1->alpha[1] += a2;
    
    gd2 = bpar->pgen2.gen_point_par;
    memcpy(gd2->alpha, bpar->prior_alpha, sizeof(gd2->alpha));
    gd2->alpha[0] += b1;
    gd2->alpha[1] += b2;

    /* 1. Find an UNCHANGED state, if it exists. */
    unsigned step = 0, a1_lo, a1_hi, max_step = MAX(a1, MAX_COUNT1 - a1);
    int16_t us = -1;
    enum fuzzy_state state;
    while (step <= max_step)
    {
        if (step <= a1)
        {
            gd1->alpha[0] = a1 + step + bpar->prior_alpha[0];
            state = binomial_quantile_est(bpar->max_points, bpar->min_dist, bpar->post_conf,
                                          bpar->beta_conf, bpar->pgen1, bpar->points1,
                                          bpar->pgen2, bpar->points2,
                                          bpar->batch_size);
            if (state == UNCHANGED)
            {
                us = a1 + step;
                break;
            }
        }
        if (step + a1 <= MAX_COUNT)
        {
            gd1->alpha[0] = a1 - step + bpar->prior_alpha[0];
            state = binomial_quantile_est(bpar->max_points, bpar->min_dist, bpar->post_conf,
                                          bpar->beta_conf, bpar->pgen1, bpar->points1,
                                          bpar->pgen2, bpar->points2,
                                          bpar->batch_size);
            if (state == UNCHANGED)
            {
                us = a1 - step;
                break;
            }
        }
        ++step;
    }
    if (us == -1) beb->unchanged[0] = beb->unchanged[1] = 0;
    else
    {
        /* found an UNCHANGED state.  Now find the bounds */
        bpar->query_state = UNCHANGED;
        beb->unchanged[0] = virtual_lower_bound(0, us, elem_is_less, bpar);
        beb->unchanged[1] = virtual_upper_bound(us, MAX_COUNT1, query_is_less, bpar);
    }
    /* Find lower bound of AMBIGUOUS region */
    bpar->query_state = AMBIGUOUS;
    beb->ambiguous[0] = virtual_lower_bound(0, beb->unchanged[0], elem_is_less, bpar);
    beb->ambiguous[1] = virtual_upper_bound(beb->unchanged[1], MAX_COUNT1, query_is_less, bpar);
}

/* test two dirichlets based on their counts. */
enum fuzzy_state cached_dirichlet_diff(unsigned *a_counts,
                                       unsigned *b_counts,
                                       struct binomial_est_params *bpar)
{
    int a1, a2, b1, b2;
    unsigned za, zb;
    find_top2_aux(a_counts, &a1, &a2, &za);
    find_top2_aux(b_counts, &b1, &b2, &zb);
    if (za >= 2 && zb >= 2
        && a_counts[a1] < MAX_COUNT1
        && a_counts[a2] < MAX_COUNT2
        && b_counts[b1] < MAX_COUNT1
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
                    initialize_est_bounds(a_counts[a2], b_counts[b1], b_counts[b2], bpar, beb);
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
        bpar->pgen1->alpha[0] = a1 + bpar->prior_alpha[0];
        bpar->pgen1->alpha[1] = a2 + bpar->prior_alpha[1];
        bpar->pgen2->alpha[0] = b1 + bpar->prior_alpha[0];
        bpar->pgen2->alpha[1] = b2 + bpar->prior_alpha[1];
        state = binomial_quantile_est(bpar->max_points, bpar->min_dist, bpar->post_conf,
                                      bpar->beta_conf, bpar->pgen1, bpar->points1,
                                      bpar->pgen2, bpar->points2,
                                      bpar->batch_size);
        return state;
    }
}
