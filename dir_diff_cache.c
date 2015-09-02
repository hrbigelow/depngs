/* provides memoization of classifying two dirichlets as same or
   different, based on their two largest components, and various
   confidence thresholds. */

#include "dir_diff_cache.h"
#include "virtual_bound.h"
#include "chunk_strategy.h"
#include "dir_points_gen.h"
#include "timer.h"

#include <float.h>
#include <assert.h>

#define MAX(a,b) ((a) < (b) ? (b) : (a))
#define MIN(a,b) ((a) < (b) ? (a) : (b))

static struct binomial_est_params g_be_par;

/* Describes estimated distance between two Dirichlets with alphas equal to:
   { A1 + p, A2 + p, p, p }
   { B1 + p, B2 + p, p, p }

   Where A1, A2, B1, B2 are integral values.  [unchanged[0],
   unchanged[1]) represents the range of values of A1 where it is
   deemed UNCHANGED.  [ambiguous[0], ambiguous[1]) are the values of
   A1 for which it is deemed AMBIGUOUS (or UNCHANGED, where this
   interval overlaps the unchanged interval).  By construction: 
   0 <= A[0] <= U[0] <= U[1] <= A[1] <= MAX_COUNT1
*/
   
/* bounds_cache[a2][b2][b1].  sets each element to have 'state' = UNSET */
void
dirichlet_diff_cache_init(struct dirichlet_diff_params dd_par,
                          struct binomial_est_params be_par,
                          struct dir_cache_params dc_par,
                          struct bam_filter_params bf_par,
                          struct bam_scanner_info *reader_buf,
                          unsigned n_max_reading,
                          unsigned long max_input_mem,
                          unsigned n_threads)
{
    g_dd_par = dd_par;
    g_be_par = be_par;

    struct dirichlet_points_gen_params pg_par = {
        .min_base_quality = bf_par.min_base_quality,
        .max_sample_points = g_be_par.max_sample_points,
        .alpha_prior = g_dd_par.prior_alpha
    };

    dirichlet_points_gen_init(pg_par);

    fprintf(stdout, "%s: Start computing confidence interval statistics.\n", timer_progress());
    // binomial_est_init(be_par, be_par.max_sample_points, n_threads);
    binomial_est_init(be_par, 12000, n_threads);
    fprintf(stdout, "%s: Finished computing confidence interval statistics.\n", timer_progress());

    dir_cache_init(dc_par);

    fprintf(stdout, "%s: Start collecting input statistics.\n", timer_progress());
    run_survey(bf_par, reader_buf, dd_par.pseudo_depth,
               dc_par.n_max_survey_loci, n_threads, n_max_reading, max_input_mem);
    fprintf(stdout, "%s: Finished collecting input statistics.\n", timer_progress());

    /* This is needed to return batch_pileup back to the beginning
       state. */
    pileup_reset_pos();
    chunk_strategy_reset();

    fprintf(stdout, "%s: Generating dirichlet point sets.\n", timer_progress());
    generate_point_sets(n_threads);
    fprintf(stdout, "%s: Finished generating dirichlet point sets.\n", timer_progress());

    fprintf(stdout, "%s: Start calculating pair-dirichlet bounds.\n", timer_progress());
    generate_est_bounds(n_threads);
    fprintf(stdout, "%s: Finished calculating pair-dirichlet bounds.\n", timer_progress());
}


void
dirichlet_diff_cache_free()
{
    binomial_est_free();
    dir_cache_free();
}


void
dir_diff_cache_thread_init()
{
    dir_points_thread_init();
}


void
dir_diff_cache_thread_free()
{
    dir_points_thread_free();
}


/* convenience function for just updating a single alpha component. */
void
set_dir_alpha_single(unsigned i, unsigned v, struct dir_points *dp)
{
    unsigned cts[4];
    memcpy(cts, dp->perm_alpha, sizeof(cts));
    cts[i] = v;
    dir_points_update_alpha(cts, NULL, dp);
}



/* mode-finding algorithm, robust to some amount of error in the
   measurement */
struct ipoint {
    unsigned x;
    double y;
    struct ipoint *left, *right, *down;
    struct binomial_est_state est;
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



/* The main distance function.  This will modify dist[0]'s alpha[0]
   parameters and thus discard its points. It will assume all 7 other
   alpha components are set as intended.  */
static struct binomial_est_state
    pair_dist_aux(struct dir_points *dp1,
                  struct dir_points *dp2,
                  unsigned a1)
{
    set_dir_alpha_single(0, a1, dp1);
    return binomial_quantile_est(dp1, dp2, GEN_POINTS_BATCH);
}
                                

/* USAGE: virtual_lower_bound(beg, end, elem_is_less, par).  Assumes
   all elements in [beg, end) are in ascending order.  Return the
   position of the left-most element that is b->query_state. */
static int 
elem_is_less(unsigned pos, void *par)
{
    struct bound_search_params *bsp = par;
    struct binomial_est_state est =
        pair_dist_aux(bsp->dist[0], bsp->dist[1], pos);
    
    return bsp->use_low_beta
        ? (est.beta_qval_lo < bsp->query_beta ? 1 : 0)
        : (est.beta_qval_hi < bsp->query_beta ? 1 : 0);
}


/* Number of highest points to explore by flanking, used in
   noisy_mode. */
#define NUMTOP 2

/* find the mode in a semi-robust way.  for the 'down' member,
   0xDEADBEEF is used to mean 'uninitialized', while NULL means 'last
   in the chain' */
static unsigned
noisy_mode(unsigned xmin, unsigned xend, void *par)
{
    /* create the two anchors */
    struct ipoint *a = malloc(sizeof(struct ipoint));
    struct ipoint *b = malloc(sizeof(struct ipoint));
    struct binomial_est_state y;
    struct bound_search_params *bsp = par;
    
    unsigned xmax = xend - 1;

    y = pair_dist_aux(bsp->dist[0], bsp->dist[1], xmin);
    a->x = xmin;
    a->y = y.beta_qval_lo;
    a->est = y;
    a->left = NULL;
    a->right = b;
    a->down = (void *)0xDEADBEEF;

    y = pair_dist_aux(bsp->dist[0], bsp->dist[1], xmax);
    b->x = xmax;
    b->y = y.beta_qval_lo;
    b->est = y;
    b->left = a;
    b->right = NULL;
    b->down = (void *)0xDEADBEEF;

    /* order the two anchors vertically */
    struct ipoint *hd = a->y > b->y ? a : b;
    struct ipoint *lo = a->y > b->y ? b : a;
    hd->down = NULL;

    INSERT_NODE_VERT(hd, lo);

    while (1) {
        unsigned i, x;
        struct ipoint *cen, *nn, *newnodes[NUMTOP * 2];
        for (i = 0; i != sizeof(newnodes) / sizeof(newnodes[0]); ++i)
            newnodes[i] = NULL;

        /* for each of the top NUMTOP nodes, create two flanking nodes. */
        for (i = 0, cen = hd; i != NUMTOP && cen; ++i) {
            if (cen->left
                && cen->left->down != (void *)0xDEADBEEF
                && (x = (cen->x + cen->left->x) / 2) != cen->x
                && x != cen->left->x) {
                nn = malloc(sizeof(struct ipoint));
                nn->x = x;
                y = pair_dist_aux(bsp->dist[0], bsp->dist[1], nn->x);
                nn->y = y.beta_qval_lo;
                nn->est = y;
                INSERT_NODE_HORZ(cen->left, nn, cen);
                nn->down = (void *)0xDEADBEEF;
                newnodes[i * 2] = nn;
            }

            if (cen->right
                && cen->right->down != (void *)0xDEADBEEF
                && (x = (cen->x + cen->right->x) / 2) != cen->x
                && x != cen->right->x) {
                nn = malloc(sizeof(struct ipoint));
                nn->x = x;
                y = pair_dist_aux(bsp->dist[0], bsp->dist[1], nn->x);
                nn->y = y.beta_qval_lo;
                nn->est = y;
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

    /* print out nodes from left to right */
    struct ipoint *tp = hd;
    while (tp->left) tp = tp->left;
    unsigned a_cts[4], b_cts[4];
    memcpy(a_cts, bsp->dist[0]->perm_alpha, sizeof(a_cts));
    memcpy(b_cts, bsp->dist[1]->perm_alpha, sizeof(b_cts));
    while (tp) {
        fprintf(stdout, 
                "%5s[%u,%u,%u,%u]\t[%u,%u,%u,%u]\t%u\t"
                "%8u\t%8u\t%20.18g\t%20.18g\n", 
                (tp == hd ? "***" : "   "),
                a_cts[0], a_cts[1], a_cts[2], a_cts[3],
                b_cts[0], b_cts[1], b_cts[2], b_cts[3],
                tp->x, 
                tp->est.n_trials, tp->est.n_success, 
                tp->est.beta_qval_lo, tp->est.beta_qval_hi);
        tp = tp->right;
    }

    /* clean up */
    struct ipoint *p;
    while (hd) {
        p = hd;
        hd = hd->down;
        free(p);
    }
    return topx;
}



/* For two dirichlet distributions A = { x+p, a2+p, p, p } and B = {
   b1+p, b2+p, p, p }, virtually find the values of x in [0, max1)
   that denote the intervals where A and B are UNCHANGED, and where
   they are AMBIGUOUS.  The expected underlying pattern as a function
   of x is some number of C, AC, A, AU, U, AU, A, AC, C.  (see enum
   fuzzy_state in binomial_est.h) */
void
initialize_est_bounds(unsigned a2, unsigned b1, unsigned b2,
                      struct bound_search_params *bsp,
                      struct binomial_est_bounds *beb)
{
    beb->ambiguous[0] = beb->ambiguous[1] = 0;
    beb->unchanged[0] = beb->unchanged[1] = 0;

    unsigned alpha1_cts[] = { 0, a2, 0, 0 };
    dir_points_update_alpha(alpha1_cts, NULL, bsp->dist[0]);

    unsigned alpha2_cts[] = { b1, b2, 0, 0 };
    dir_points_update_alpha(alpha2_cts, NULL, bsp->dist[1]);

    /* Find Mode.  (consider [0, xmode) and [xmode, cache.max1) as the
       upward and downward phase intervals */
    unsigned xmode = noisy_mode(0, g_dd_par.pseudo_depth, bsp);

    bsp->use_low_beta = 0;
    bsp->query_beta = 1.0 - g_be_par.post_confidence;
    beb->ambiguous[0] = virtual_lower_bound(0, xmode, elem_is_less, bsp);

    /* find position in the virtual array of distances [xmode, g_dd_par.pseudo_depth) */
    /* the range [xmode, g_dd_par.pseudo_depth) is monotonically
       DECREASING. virtual_upper_bound assumes monotonically
       increasing range, so we must use elem_is_less as the function
       rather than query_is_less. */
    beb->ambiguous[1] = virtual_upper_bound(xmode, g_dd_par.pseudo_depth, elem_is_less, bsp);

    bsp->use_low_beta = 1;
    bsp->query_beta = g_be_par.post_confidence;
    beb->unchanged[0] = virtual_lower_bound(beb->ambiguous[0], xmode, elem_is_less, bsp);

    beb->unchanged[1] = virtual_upper_bound(xmode, beb->ambiguous[1], elem_is_less, bsp);
    
}


/* attempt to retrieve the fuzzy_state alpha1 vs alpha2 distance from
   the dir_cache, or calculate it if not found */
enum fuzzy_state
cached_dirichlet_diff(unsigned *alpha1, unsigned *alpha2,
                      struct bound_search_params *bsp,
                      unsigned *cache_hit)
{
    unsigned perm[4], *perm_used;
    enum fuzzy_state state = UNCHANGED;
    enum diff_cache_status status;
    
    status = dir_cache_try_get_diff(alpha1, alpha2, perm, &state);
    if (status == CACHE_HIT) {
        *cache_hit = 1;
        return state;
    } else {
        perm_used = (status == CACHE_MISS_PERM) ? perm : NULL;
        dir_points_update_alpha(alpha1, perm_used, bsp->dist[0]);
        dir_points_update_alpha(alpha2, perm_used, bsp->dist[1]);
        struct binomial_est_state est =
            binomial_quantile_est(bsp->dist[0], bsp->dist[1], GEN_POINTS_BATCH);
        state = est.state;
        return state;
    }
}
