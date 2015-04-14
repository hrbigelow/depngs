/* provides memoization of classifying two dirichlets as same or
   different, based on their two largest components, and various
   confidence thresholds. */

#include "binomial_est.h"
#include "virtual_bound.h"
#include "dirichlet_diff_cache.h"
#include "dirichlet_points_gen.h"

#include <math.h>
#include <string.h>
#include <pthread.h>
#include <float.h>
#include <inttypes.h>

#define MAX_COUNT1 50
#define MAX_COUNT2 10

#define MAX(a,b) ((a) < (b) ? (b) : (a))
#define MIN(a,b) ((a) < (b) ? (a) : (b))

pthread_mutex_t set_flag_mtx;
struct binomial_est_bounds ***bounds_cache, **bounds_cache_buf2, *bounds_cache_buf;

/* bounds_cache[a2][b2][b1] */
void dirichlet_diff_init(unsigned max1, unsigned max2)
{
    pthread_mutex_init(&set_flag_mtx, NULL);
    bounds_cache_buf = malloc(max2 * max2 * max1 * sizeof(struct binomial_est_bounds));
    bounds_cache_buf2 = malloc(max2 * max2 * sizeof(struct binomial_est_bounds *));
    bounds_cache = malloc(max2 * sizeof(struct binomial_est_bounds **));

    unsigned i, j, k, l;
    unsigned M = max2 * max2;
    for (i = 0, j = 0; i != M; ++i, j += max1)
        bounds_cache_buf2[i] = bounds_cache_buf + j;

    for (k = 0, l = 0; k != max2; ++k, l += max2)
        bounds_cache[k] = bounds_cache_buf2 + l;
}


void dirichlet_diff_free()
{
    free(bounds_cache);
    free(bounds_cache_buf);
    free(bounds_cache_buf2);
    pthread_mutex_destroy(&set_flag_mtx);
}


void alloc_distrib_points(struct distrib_points *dpts,
                          unsigned max_sample_points)
{
    unsigned msp = max_sample_points;
    dpts->pgen = (struct points_gen){ 
        malloc(sizeof(struct dir_points_par)),
        gen_dirichlet_points_wrapper, 
        malloc(sizeof(struct calc_post_to_dir_par)),
        calc_post_to_dir_ratio
    };
    ((struct dir_points_par *)dpts->pgen.point_par)->randgen = 
        gsl_rng_alloc(gsl_rng_taus);
    dpts->points = (struct points_buf){ (POINT *)malloc(sizeof(POINT) * msp), 0, msp };
    dpts->weights = (struct weights_buf){ (double *)malloc(sizeof(double) * msp), 0, msp };
}


void free_distrib_points(struct distrib_points *dpts)
{
    gsl_rng_free(((struct dir_points_par *)dpts->pgen.point_par)->randgen);
    free((struct dir_points_par *)dpts->pgen.point_par);
    free((struct calc_post_to_dir_par *)dpts->pgen.weight_par);
    free(dpts->points.buf);
    free(dpts->weights.buf);
}


/* update dpts parameters and discard existing points.  if alpha is
   identical to existing parameters, do nothing. */
void set_dirichlet_alpha(struct distrib_points *dpts, double *alpha)
{
    struct dir_points_par *dp = dpts->pgen.point_par;
    unsigned i, change = 0;
    for (i = 0; i != NUM_NUCS; ++i)
        if (dp->alpha[i] != alpha[i]) change = 1;

    if (change)
    {
        memcpy(dp->alpha, alpha, sizeof(dp->alpha));
        dpts->points.size = 0;
    }
}

/* set a single alpha, and discard points if there is a change */
void set_dirichlet_alpha_single(struct distrib_points *dpts, unsigned i, double v)
{
    struct dir_points_par *dp = dpts->pgen.point_par;
    if (dp->alpha[i] != v)
    {
        dp->alpha[i] = v;
        dpts->points.size = 0;
    }
}

/* The main distance function.  This will modify dist[0]'s alpha
   parameters and thus discard its points. */
struct binomial_est_state pair_dist_aux(unsigned a1, void *par)
{
    struct binomial_est_params *b = par;
    struct posterior_settings *ps = b->pset;

    set_dirichlet_alpha_single(b->dist[0], 0, a1 + ps->prior_alpha[0]);

    struct binomial_est_state est;
    est = binomial_quantile_est(ps->max_sample_points, 
                                ps->min_dist,
                                ps->post_confidence, 
                                ps->beta_confidence,
                                b->dist[0]->pgen,
                                &b->dist[0]->points, 
                                b->dist[1]->pgen,
                                &b->dist[1]->points,
                                b->batch_size);
    return est;
}


void read_diststats_line(FILE *fh, 
                         unsigned *a2,
                         unsigned *b1,
                         unsigned *b2,
                         struct binomial_est_bounds *beb)
{
    int n;
    n = fscanf(fh, 
               "%u\t%u\t%u\t%"SCNd16"\t%"SCNd16"\t%"SCNd16"\t%"SCNd16"\n", 
               a2, b1, b2,
               &beb->ambiguous[0],
               &beb->unchanged[0],
               &beb->unchanged[1],
               &beb->ambiguous[1]);
    if (n != 7)
    {
        fprintf(stderr, "error reading diststats line at %s: %u\n", __FILE__, __LINE__);
        exit(1);
    }
}


void write_diststats_line(FILE *fh,
                          unsigned a2,
                          unsigned b1,
                          unsigned b2,
                          struct binomial_est_bounds *beb)
{
    fprintf(fh,
            "%i\t%i\t%i\t%i\t%i\t%i\t%i\n", 
            a2, b1, b2,
            beb->ambiguous[0], beb->unchanged[0], 
            beb->unchanged[1], beb->ambiguous[1]);
}

/* parse diststats header, initializing fields of pset and max1 and max2  */
void parse_diststats_header(FILE *diststats_fh, 
                            struct posterior_settings *pset,
                            unsigned *max1, 
                            unsigned *max2)
{
    int n;
    float prior_alpha;
    n = fscanf(diststats_fh, 
               "# max1: %u\n"
               "# max2: %u\n"
               "# max_sample_points: %zu\n"
               "# min_dist: %lf\n"
               "# post_confidence: %lf\n"
               "# beta_confidence: %lf\n"
               "# prior_alpha: %f\n",
               max1, max2, 
               &pset->max_sample_points,
               &pset->min_dist,
               &pset->post_confidence,
               &pset->beta_confidence,
               &prior_alpha);
    if (n != 7)
    {
        fprintf(stderr, "error parsing diststats header at %s: %ul\n", __FILE__, __LINE__);
        exit(1);
    }
    unsigned i;
    for (i = 0; i != NUM_NUCS; ++i)
        pset->prior_alpha[i] = prior_alpha;
}

void write_diststats_header(FILE *diststats_fh,
                            struct posterior_settings pset,
                            unsigned max1,
                            unsigned max2)
{
    fprintf(diststats_fh, 
            "# max1: %u\n"
            "# max2: %u\n"
            "# max_sample_points: %zu\n"
            "# min_dist: %lf\n"
            "# post_confidence: %lf\n"
            "# beta_confidence: %lf\n"
            "# prior_alpha: %f\n",
            max1, max2, 
            pset.max_sample_points,
            pset.min_dist,
            pset.post_confidence,
            pset.beta_confidence,
            pset.prior_alpha[0]);
}


/* initialize internal bounds_cache */
void parse_diststats_body(FILE *diststats_fh, unsigned max1, unsigned max2)
{
    struct binomial_est_bounds beb;
    unsigned a2, b1, b2;
    while (! feof(diststats_fh))
    {
        read_diststats_line(diststats_fh, &a2, &b1, &b2, &beb);
        if (a2 < max2 && b2 < max2 && b1 < max1)
            bounds_cache[a2][b2][b1] = beb;
        else
        {
            fprintf(stderr, 
                    "Error at %s: %u\n"
                    "Indices exceed limits.  A2 = %u, B2 = %u, B1 = %u, MAX1 = %u, MAX2 = %u\n",
                    __FILE__, __LINE__,
                    a2, b2, b1, max1, max2);
            exit(1);
        }
    }
}


/* mode-finding algorithm, robust to some amount of error in the
   measurement */
struct ipoint {
    unsigned x;
    double y;
    struct ipoint *left, *right, *down;
};


#define INSERT_NODE_HORZ(L, M, R)  \
    do {                           \
        struct ipoint              \
            *Lp = (L),             \
            *Mp = (M),             \
            *Rp = (R);             \
        Mp->left = Lp;             \
        Mp->right = Rp;            \
        Lp->right = Mp;            \
        Rp->left = Mp;             \
    } while (0)                    \
        

/* Preconditions:  
   1. HEAD != NULL
   2. HEAD->down == NULL || HEAD->y > HEAD->down->y
   3. N != NULL and N is not in the chain
   4. N->x does not appear as an x value in any of the nodes in the chain.

   Postcondition:
   5. N appears in the chain, and the chain is well-ordered vertically.
   6. HEAD points to the true head of the updated list, which may be N.
*/
#define INSERT_NODE_VERT(HEAD, N)                               \
    do {                                                        \
        struct ipoint pre = { 0, DBL_MAX, NULL, NULL, (HEAD) }; \
        struct ipoint *p = &pre;                                \
        while (p && (p->down ? p->down->y : 0) > (N)->y)        \
            p = p->down;                                        \
        (N)->down = p->down;                                    \
        p->down = (N);                                          \
        (HEAD) = pre.down;                                      \
    } while (0)                                                 \
        



/* Number of highest points to explore by flanking, used in
   noisy_mode. */
#define NUMTOP 2

/* find the mode in a semi-robust way.  for the 'down' member,
   0xDEADBEEF is used to mean 'uninitialized', while NULL means 'last
   in the chain'
 */
unsigned noisy_mode(unsigned xmin, unsigned xend, void *bpar)
{
    /* create the two anchors */
    struct ipoint *a = malloc(sizeof(struct ipoint));
    struct ipoint *b = malloc(sizeof(struct ipoint));
    struct binomial_est_state y;
    unsigned xmax = xend - 1;

    y = pair_dist_aux(xmin, bpar);
    a->x = xmin, a->y = y.beta_qval_lo, a->left = NULL, a->right = b, a->down = (void *)0xDEADBEEF;

    y = pair_dist_aux(xmax, bpar);
    b->x = xmax, b->y = y.beta_qval_lo, b->left = a, b->right = NULL, b->down = (void *)0xDEADBEEF;

    /* order the two anchors vertically */
    struct ipoint *hd = a->y > b->y ? a : b;
    struct ipoint *lo = a->y > b->y ? b : a;
    hd->down = NULL;

    INSERT_NODE_VERT(hd, lo);

    while (1)
    {
        unsigned i, x;
        struct ipoint *cen, *nn, *newnodes[NUMTOP * 2];
        for (i = 0; i != sizeof(newnodes) / sizeof(newnodes[0]); ++i)
            newnodes[i] = NULL;

        /* for each of the top NUMTOP nodes, create two flanking nodes. */
        for (i = 0, cen = hd; i != NUMTOP && cen; ++i)
        {
            if (cen->left
                && cen->left->down != (void *)0xDEADBEEF
                && (x = (cen->x + cen->left->x) / 2) != cen->x
                && x != cen->left->x)
            {
                nn = malloc(sizeof(struct ipoint));
                nn->x = x;
                y = pair_dist_aux(nn->x, bpar);
                nn->y = y.beta_qval_lo;
                INSERT_NODE_HORZ(cen->left, nn, cen);
                nn->down = (void *)0xDEADBEEF;
                newnodes[i * 2] = nn;
            }

            if (cen->right
                && cen->right->down != (void *)0xDEADBEEF
                && (x = (cen->x + cen->right->x) / 2) != cen->x
                && x != cen->right->x)
            {
                nn = malloc(sizeof(struct ipoint));
                nn->x = x;
                y = pair_dist_aux(nn->x, bpar);
                nn->y = y.beta_qval_lo;
                INSERT_NODE_HORZ(cen, nn, cen->right);
                nn->down = (void *)0xDEADBEEF;
                newnodes[i * 2 + 1] = nn;
            }

            cen = cen->down;
        }
        /* Now, up to 4 new nodes are created.  They will be the
           immediate neighbors of the top two nodes */
        for (i = 0; i != NUMTOP * 2; ++i)
            if (newnodes[i]) INSERT_NODE_VERT(hd, newnodes[i]);

        /* We want these to represent distance to neighbor along x
           axis. Where there is no neighbor, consider it to be '1' */
        unsigned lx = (int)hd->x - (hd->left ? hd->left->x : -1);
        unsigned rx = (hd->right ? hd->right->x : xend) - hd->x;

        /* stop iff both neighbors are as close as possible */
        if (lx == 1 && rx == 1) break;
    }
    unsigned topx = hd->x;

    /* struct binomial_est_params *bb = bpar; */
    /* struct dir_points_par  */
    /*     *d0 = bb->dist[0]->pgen.point_par, */
    /*     *d1 = bb->dist[1]->pgen.point_par; */

    /* fprintf(stderr, "MODE: %5.3g\t%5.3g\t%5.3g\t%5.3g\t%i\t%7.4g\n",  */
    /*         d0->alpha[0], d0->alpha[1], d1->alpha[0], d1->alpha[1], */
    /*         hd->x, hd->y); */

    /* clean up */
    struct ipoint *p;
    while (hd)
    {
        p = hd;
        hd = hd->down;
        free(p);
    }
    return topx;
}


/* bounds_cache[a2][b2][b1] = a description of the matrix of distance
   categories.  */
// struct binomial_est_bounds bounds_cache[MAX_COUNT2][MAX_COUNT2][MAX_COUNT1];



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


/* Given a[4], b[4], and lim[4], find a permutation of {(a[p_1],
   b[p_1]), (a[p_2], b[p_2]), (a[p_3], b[p_e]), (a[p_4], b[p_4]) }
   such that both a[p_i] and b[p_i] are less than lim[p_i], if such
   permutation exists. 

   Modification: Given lim[2] == 1 and lim[3] == 1, The permutation of
   the last two elements does not matter.  So, we only need to specify
   p_1 and p_2 in this case.

   Also, if there is no permutation that suffices, set a flag
   appropriately.

   Assumes that lim is descending.  lim[i] >= lim[i+1]
*/
void find_cacheable_permutation(const unsigned *a, const unsigned *b, 
                                const unsigned *lim,
                                unsigned char *permutation, 
                                unsigned *perm_found)
{
    *perm_found = 1;
    /* mpi[i] (max permutation index) is the maximum position in the
       permutation (described above) that the (a, b) pair may attain
       and still stay below the limits. */
    int mpi[] = { -1, -1, -1, -1 },
        i, j,
        mv = 0, mi = 0; /* min value, min_index */

    for (i = 0; i != 4; ++i)
        for (j = 0; j != 4; ++j)
        {
            if (a[i] >= lim[j] || b[i] >= lim[j]) break;
            mpi[i] = j;
        }

    /* p[i], (the permutation of mpi), s.t. mpi[p[i]] <= mpi[p[i+1]] */
    unsigned p[4];

    /* flag[j] = 1 means the value j has been placed in the permutation */
    unsigned char flag[] = { 0, 0, 0, 0 };

    for (i = 0; i != 4; ++i)
    {
        mv = 5;
        for (j = 0; j != 4; ++j)
            if (! flag[j] && mpi[j] < mv)
                mv = mpi[j], mi = j;
        p[i] = mi;
        flag[mi] = 1;
    }
    for (i = 0; i != 4; ++i)
        if (mpi[p[i]] < i) { *perm_found = 0; break; }

    permutation[0] = p[0];
    permutation[1] = p[1];
}

                                

/* find top 2 components in el, store their indices in i1 and i2.
   Count the number of occurrences of min_val in z */
void find_top2_aux(unsigned *el, int *i1, int *i2, unsigned *z)
{
    unsigned m1 = 0, m2 = 0;
    *i1 = 0, *i2 = 0;
    *z = 0;
    unsigned i;
    for (i = 0; i != NUM_NUCS; ++i)
    {
        if (m1 < el[i]) *i1 = i, m1 = el[i];
        if (el[i] == 0) ++*z;
    }
    /* set i2 to an arbitrary value that is unequal to i1 */
    *i2 = (NUM_NUCS - 1) - *i1;
    m2 = el[*i2];
    for (i = 0; i != NUM_NUCS; ++i)
        if (i != *i1 && m2 < el[i])
            *i2 = i, m2 = el[i];
}



/* depends on the ordering of enum fuzzy_state in a gradient from
   CHANGED to UNCHANGED.
   */
int query_is_less(unsigned pos, void *par)
{
    struct binomial_est_params *b = par;
    struct binomial_est_state est = pair_dist_aux(pos, par);
    return b->use_low_beta
        ? (b->query_beta < est.beta_qval_lo ? 1 : 0)
        : (b->query_beta < est.beta_qval_hi ? 1 : 0);
}

/* USAGE: virtual_lower_bound(beg, end, elem_is_less, par).  Assumes
   all elements in [beg, end) are in ascending order.  Return the
   position of the left-most element that is b->query_state. */
int elem_is_less(unsigned pos, void *par)
{
    struct binomial_est_params *b = par;
    struct binomial_est_state est = pair_dist_aux(pos, par);
    return b->use_low_beta
        ? (est.beta_qval_lo < b->query_beta ? 1 : 0)
        : (est.beta_qval_hi < b->query_beta ? 1 : 0);
}


#define CACHE_GEN_POINTS_BATCH 20




/* For two dirichlet distributions A = { x+p, a2+p, p, p } and B = {
   b1+p, b2+p, p, p }, virtually find the values of x in [0,
   MAX_COUNT1) that denote the intervals where A and B are UNCHANGED,
   and where they are AMBIGUOUS.  The expected underlying pattern as a
   function of x is some number of C, AC, A, AU, U, AU, A, AC, C.
   (see enum fuzzy_state in binomial_est.h) */
void initialize_est_bounds(unsigned a2, unsigned b1, unsigned b2,
                           struct binomial_est_params *bpar,
                           struct binomial_est_bounds *beb)
{
    beb->ambiguous[0] = beb->ambiguous[1] = 0;
    beb->unchanged[0] = beb->unchanged[1] = 0;

    struct posterior_settings *ps = bpar->pset;

    double alpha[NUM_NUCS];
    memcpy(alpha, bpar->pset->prior_alpha, sizeof(alpha));
    alpha[1] += a2;
    set_dirichlet_alpha(bpar->dist[0], alpha);

    memcpy(alpha, bpar->pset->prior_alpha, sizeof(alpha));
    alpha[0] += b1;
    alpha[1] += b2;
    set_dirichlet_alpha(bpar->dist[1], alpha);

    /* Find Mode.  (consider [0, xmode) and [xmode, bpar->max1) as the
       upward and downward phase intervals */
    unsigned xmode = noisy_mode(0, bpar->max1, bpar);

    bpar->use_low_beta = 0;
    bpar->query_beta = 1.0 - ps->post_confidence;
    beb->ambiguous[0] = virtual_lower_bound(0, xmode, elem_is_less, bpar);
    beb->ambiguous[1] = virtual_upper_bound(xmode, bpar->max1, elem_is_less, bpar);

    bpar->use_low_beta = 1;
    bpar->query_beta = ps->post_confidence;
    beb->unchanged[0] = virtual_lower_bound(beb->ambiguous[0], xmode, elem_is_less, bpar);
    beb->unchanged[1] = virtual_upper_bound(xmode, beb->ambiguous[1], elem_is_less, bpar);
    
}

#if 0
/* Does the same as print_bounds, but calculates the bounds using
   binary search. */
void print_beb_bounds(struct binomial_est_params *bpar)
{
    unsigned a1, a2, b1, b2;
    char states[] = "CRAFU";
    struct binomial_est_bounds beb;
    // struct binomial_est_state est;
    // struct posterior_settings *ps = bpar->pset;
    // double alpha[NUM_NUCS];

    /* a1 >= a2, b1 >= b2*/
    for (b1 = 0; b1 != bpar->max1; ++b1)
    // for (b1 = 30; b1 != 40; ++b1)
    {
        for (a2 = 0; a2 != bpar->max2 && a2 <= a1; ++a2)
            // for (a2 = 0; a2 != 1; ++a2)
        {
            for (b2 = 0; b2 != bpar->max2 && b2 <= b1; ++b2)
                // for (b2 = 0; b2 != 5 && b2 <= b1; ++b2)
            {
                initialize_est_bounds(a2, b1, b2, bpar, &beb);
                for (a1 = a2; a1 != bpar->max1; ++a1)
                    fprintf(stdout, "PBB\t%i\t%i\t%i\t%i\t0.0\t0.0\t%c\n", a1, a2, b1, b2,
                            a1 < beb.ambiguous[0] ? states[CHANGED]
                            : (a1 < beb.unchanged[0] ? states[AMBIGUOUS]
                               : (a1 < beb.unchanged[1] ? states[UNCHANGED]
                                  : (a1 < beb.ambiguous[1] ? states[AMBIGUOUS]
                                     : states[CHANGED]))));
#if 0                
                for (a1 = a2; a1 != bpar->max1; ++a1)
                {
                    memcpy(alpha, bpar->pset->prior_alpha, sizeof(alpha));
                    alpha[0] += a1;
                    alpha[1] += a2;
                    set_dirichlet_alpha(bpar->dist[0], alpha);
                
                    memcpy(alpha, bpar->pset->prior_alpha, sizeof(alpha));
                    alpha[0] += b1;
                    alpha[1] += b2;
                    set_dirichlet_alpha(bpar->dist[1], alpha);
                
                    est = binomial_quantile_est(ps->max_sample_points, 
                                                ps->min_dist,
                                                ps->post_confidence, 
                                                ps->beta_confidence,
                                                bpar->dist[0]->pgen,
                                                &bpar->dist[0]->points, 
                                                bpar->dist[1]->pgen,
                                                &bpar->dist[1]->points,
                                                bpar->batch_size);
                    fprintf(stdout, "PB\t%i\t%i\t%i\t%i\t%7.5g\t%7.5g\t%c\n", a1, a2, b1, b2, 
                            est.beta_qval_lo, est.beta_qval_hi, states[est.state]);
                }
#endif
            }
        }
    }
    fflush(stdout);
}
#endif

#if 0
void print_bounds(struct binomial_est_params *bpar)
{
    unsigned a1, a2, b1, b2;
    struct binomial_est_state est;
    struct posterior_settings *ps = bpar->pset;
    char states[] = "CRAFU";
    double alpha[NUM_NUCS];
    /* a1 >= a2, b1 >= b2*/
    
    for (a1 = 0; a1 != bpar->max1; ++a1)
    {
        for (b1 = 0; b1 != bpar->max1; ++b1)
        // for (b1 = 17; b1 != 18; ++b1)
        {
            for (a2 = 0; a2 != 2 && a2 <= a1; ++a2)
                // for (a2 = 8; a2 != 9; ++a2)
            {
                for (b2 = 0; b2 != bpar->max2 && b2 <= b1; ++b2)
                    // for (b2 = 8; b2 != 9; ++b2)
                {
                    memcpy(alpha, bpar->pset->prior_alpha, sizeof(alpha));
                    alpha[0] += a1;
                    alpha[1] += a2;
                    set_dirichlet_alpha(bpar->dist[0], alpha);

                    memcpy(alpha, bpar->pset->prior_alpha, sizeof(alpha));
                    alpha[0] += b1;
                    alpha[1] += b2;
                    set_dirichlet_alpha(bpar->dist[1], alpha);

                    est = binomial_quantile_est(ps->max_sample_points, 
                                                ps->min_dist,
                                                ps->post_confidence, 
                                                ps->beta_confidence,
                                                bpar->dist[0]->pgen,
                                                &bpar->dist[0]->points, 
                                                bpar->dist[1]->pgen,
                                                &bpar->dist[1]->points,
                                                bpar->batch_size);
                    fprintf(stdout, "PB\t%i\t%i\t%i\t%i\t%7.5g\t%7.5g\t%c\n", a1, a2, b1, b2, 
                            est.beta_qval_lo, est.beta_qval_hi, states[est.state]);
                }
            }
        }
    }
    fflush(stdout);
}
#endif


/* test two dirichlets based on their counts. */
enum fuzzy_state cached_dirichlet_diff(unsigned *a_counts,
                                       unsigned *b_counts,
                                       struct binomial_est_params *bpar,
                                       unsigned *cache_hit)
{
    unsigned lim[] = { bpar->max1, bpar->max2, 1, 1 };
    unsigned char perm[2];
    
    find_cacheable_permutation(a_counts, b_counts, lim, perm, cache_hit);
    
    enum fuzzy_state state;
    if (*cache_hit)
    {
        /* fprintf(stderr, "C***: %u,%u,%u,%u\t%u,%u,%u,%u\n", */
        /*         a_counts[0], a_counts[1], a_counts[2], a_counts[3], */
        /*         b_counts[0], b_counts[1], b_counts[2], b_counts[3]); */
        
        unsigned i1 = perm[0], i2 = perm[1];

        struct binomial_est_bounds *beb = 
            &bounds_cache[a_counts[i2]][b_counts[i2]][b_counts[i1]];
        
        state = locate_cell(beb, a_counts[i1]);
    }
    else
    {
        /* fprintf(stderr, "N---: %u,%u,%u,%u\t%u,%u,%u,%u\n", */
        /*         a_counts[0], a_counts[1], a_counts[2], a_counts[3], */
        /*         b_counts[0], b_counts[1], b_counts[2], b_counts[3]); */
        
        /* need to test the bounds organically, with no caching */
        struct posterior_settings *ps = bpar->pset;
        double alpha[NUM_NUCS];
        memcpy(alpha, bpar->pset->prior_alpha, sizeof(alpha));
        unsigned i;
        for (i = 0; i != NUM_NUCS; ++i) alpha[i] += a_counts[i];
        set_dirichlet_alpha(bpar->dist[0], alpha);
        
        memcpy(alpha, bpar->pset->prior_alpha, sizeof(alpha));
        for (i = 0; i != NUM_NUCS; ++i) alpha[i] += b_counts[i];
        set_dirichlet_alpha(bpar->dist[1], alpha);

        struct binomial_est_state est = 
            binomial_quantile_est(ps->max_sample_points, 
                                  ps->min_dist, 
                                  ps->post_confidence,
                                  ps->beta_confidence, 
                                  bpar->dist[0]->pgen,
                                  &bpar->dist[0]->points, 
                                  bpar->dist[1]->pgen,
                                  &bpar->dist[1]->points,
                                  bpar->batch_size);
        state = est.state;
    }
    return state;
}
