/* provides memoization of classifying two dirichlets as same or
   different, based on their two largest components, and various
   confidence thresholds. */

#include "dir_diff_cache.h"
#include "virtual_bound.h"
#include "chunk_strategy.h"
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
   






/* use to find a joint permutation between a and b, such that a2, b1,
   b2 could be packable. */
static unsigned bounds_perm_limits[] = { 
    (1<<24) + 1, (1<<20) + 1, 1, 1 
};
    



struct alpha_packed {
    unsigned a0 :12; /* 4,096 */
    unsigned a1 :10; /* 1,024 */
    unsigned a2 :6;  /*    64 */
    unsigned a3 :4;  /*    16 */
};


union alpha_key {
    struct alpha_packed c[2];
    uint64_t raw;
};


static unsigned primary_cache_size;

void print_primary_cache_size()
{
    fprintf(stderr, "Primary cache size: %u\n", primary_cache_size);
}

#define INIT_ALPHA_KEY(a1, a2, prm)                                         \
    { .c = {                                                            \
        { (a1)[(prm)[0]], (a1)[(prm)[1]], (a1)[(prm)[2]], (a1)[(prm)[3]] }, \
        { (a2)[(prm)[0]], (a2)[(prm)[1]], (a2)[(prm)[2]], (a2)[(prm)[3]] } \
        }                                                               \
    }                                                                   \
        

/* number of loci to process before freezing the bounds hash. */
#define N_LOCI_TO_FREEZE 1e8

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

/* update dpts parameters and discard existing points.  if alpha is
   identical to existing parameters, do nothing. */
void
update_points_gen_params(struct distrib_points *dpts,
                         unsigned *alpha_counts,
                         unsigned *perm)
{
    struct points_gen_par *dp = dpts->pgen.points_gen_par;
    unsigned i, change = 0;
    for (i = 0; i != NUM_NUCS; ++i) {
        if (dp->alpha_counts[i] != alpha_counts[perm[i]]
            || dp->alpha_perm[i] != perm[i]) change = 1;
        dp->alpha_counts[i] = alpha_counts[perm[i]];
        dp->alpha_perm[i] = perm[i];
    }

    if (change) {
        dpts->points.size = 0;
        dpts->weights.size = 0;
    }
}

/* set a single alpha, and discard points if there is a change */
void set_dirichlet_alpha_single(struct distrib_points *dpts, unsigned i, unsigned v)
{
    struct points_gen_par *dp = dpts->pgen.points_gen_par;
    if (dp->alpha_counts[i] != v) {
        dp->alpha_counts[i] = v;
        dpts->points.size = 0;
        dpts->weights.size = 0;
    }
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
        







/* estimate distance between two initialized 4D dirichlet
   distributions. Assume bpar->dist[0] and bpar->dist[1] have already
   had alpha_counts initialized, and have their size field set
   appropriately (to zero if the points buffer is invalid). Use
   caching of individual dirichlet distributions whenever possible. */
struct binomial_est_state 
get_est_state(struct bound_search_params *bpar)
{
    khiter_t itr;
    unsigned i;
    for (i = 0; i != 2; ++i) {
        khint64_t key;
        struct points_gen_par *pgp = bpar->dist[i]->pgen.points_gen_par;
        unsigned *cts = pgp->alpha_counts;
        unsigned packable = init_alpha_packed_large(cts, &key);
        struct points_buf *pb = &bpar->dist[i]->points;
        if (pb->size == 0 
            && packable
            && (itr = kh_get(points_h, g_points_hash, key)) != kh_end(g_points_hash)
            && kh_exist(g_points_hash, itr)) {
            pb->p = kh_val(g_points_hash, itr);
            pb->size = g_be_par.max_sample_points;
        } else {
            pb->p = pb->buf;
            pb->size = 0;
        }
    }

    struct binomial_est_state rval =
        binomial_quantile_est(bpar->dist[0]->pgen,
                              &bpar->dist[0]->points, 
                              bpar->dist[1]->pgen,
                              &bpar->dist[1]->points);
    return rval;
}


/* The main distance function.  This will modify dist[0]'s alpha[0]
   parameters and thus discard its points. It will assume all 7 other
   alpha components are set as intended.  */
static struct binomial_est_state
pair_dist_aux(unsigned a1, void *par)
{
    struct bound_search_params *bpar = par;
    set_dirichlet_alpha_single(bpar->dist[0], 0, a1);
    return get_est_state(bpar);
}
                                

/* depends on the ordering of enum fuzzy_state in a gradient from
   CHANGED to UNCHANGED.
   */
/*
static int
query_is_less(unsigned pos, void *par)
{
    struct bound_search_params *b = par;
    struct binomial_est_state est = pair_dist_aux(pos, par);
    return b->use_low_beta
        ? (b->query_beta < est.beta_qval_lo ? 1 : 0)
        : (b->query_beta < est.beta_qval_hi ? 1 : 0);
}
*/

/* USAGE: virtual_lower_bound(beg, end, elem_is_less, par).  Assumes
   all elements in [beg, end) are in ascending order.  Return the
   position of the left-most element that is b->query_state. */
static int 
elem_is_less(unsigned pos, void *par)
{
    struct bound_search_params *b = par;
    struct binomial_est_state est = pair_dist_aux(pos, par);
    return b->use_low_beta
        ? (est.beta_qval_lo < b->query_beta ? 1 : 0)
        : (est.beta_qval_hi < b->query_beta ? 1 : 0);
}


/* Number of highest points to explore by flanking, used in
   noisy_mode. */
#define NUMTOP 2

/* find the mode in a semi-robust way.  for the 'down' member,
   0xDEADBEEF is used to mean 'uninitialized', while NULL means 'last
   in the chain'
*/
static unsigned
noisy_mode(unsigned xmin, unsigned xend, void *bpar)
{
    /* create the two anchors */
    struct ipoint *a = malloc(sizeof(struct ipoint));
    struct ipoint *b = malloc(sizeof(struct ipoint));
    struct binomial_est_state y;
    unsigned xmax = xend - 1;

    y = pair_dist_aux(xmin, bpar);
    a->x = xmin;
    a->y = y.beta_qval_lo;
    a->est = y;
    a->left = NULL;
    a->right = b;
    a->down = (void *)0xDEADBEEF;

    y = pair_dist_aux(xmax, bpar);
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
                y = pair_dist_aux(nn->x, bpar);
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
                y = pair_dist_aux(nn->x, bpar);
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
    struct bound_search_params *bsp = bpar;
    memcpy(a_cts, ((struct points_gen_par *)bsp->dist[0]->pgen.points_gen_par)->alpha_counts, sizeof(a_cts));
    memcpy(b_cts, ((struct points_gen_par *)bsp->dist[1]->pgen.points_gen_par)->alpha_counts, sizeof(b_cts));
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
                      struct bound_search_params *bpar,
                      struct binomial_est_bounds *beb)
{
    beb->ambiguous[0] = beb->ambiguous[1] = 0;
    beb->unchanged[0] = beb->unchanged[1] = 0;

    unsigned perm[] = { 0, 1, 2, 3 };
    unsigned alpha1_cts[] = { 0, a2, 0, 0 };
    update_points_gen_params(bpar->dist[0], alpha1_cts, perm);

    unsigned alpha2_cts[] = { b1, b2, 0, 0 };
    update_points_gen_params(bpar->dist[1], alpha2_cts, perm);

    /* Find Mode.  (consider [0, xmode) and [xmode, cache.max1) as the
       upward and downward phase intervals */
    unsigned xmode = noisy_mode(0, g_dd_par.pseudo_depth, bpar);

    bpar->use_low_beta = 0;
    bpar->query_beta = 1.0 - g_be_par.post_confidence;
    beb->ambiguous[0] = virtual_lower_bound(0, xmode, elem_is_less, bpar);

    /* find position in the virtual array of distances [xmode, g_dd_par.pseudo_depth) */
    /* the range [xmode, g_dd_par.pseudo_depth) is monotonically
       DECREASING. virtual_upper_bound assumes monotonically
       increasing range, so we must use elem_is_less as the function
       rather than query_is_less. */
    beb->ambiguous[1] = virtual_upper_bound(xmode, g_dd_par.pseudo_depth, elem_is_less, bpar);

    bpar->use_low_beta = 1;
    bpar->query_beta = g_be_par.post_confidence;
    beb->unchanged[0] = virtual_lower_bound(beb->ambiguous[0], xmode, elem_is_less, bpar);

    beb->unchanged[1] = virtual_upper_bound(xmode, beb->ambiguous[1], elem_is_less, bpar);
    
}


/* attempt to retrieve a cached distance classification of { a1, a2,
   0, 0 } and { b1, b2, 0, 0 }.  if the cache contained this entry,
   set *state and return 1.  otherwise, return 0 */
int
get_fuzzy_state(unsigned a1, unsigned a2, 
                unsigned b1, unsigned b2,
                enum fuzzy_state *state)
{
    struct binomial_est_bounds *beb = dir_cache_try_get_bounds(a2, b1, b2);
    if (beb) {
        *state = locate_cell(beb, a1);
        return 1;
    } else 
        return 0;
}


/* test two dirichlets based on their counts. if the pattern of counts
   is cacheable, 'cacheable' is set to 1, and the cache is queried for
   the entry, and cache_was_set is set to 1 if found. */
enum fuzzy_state
cached_dirichlet_diff(unsigned *a, unsigned *b,
                      struct bound_search_params *bpar,
                      unsigned *cache_was_set)
{
    unsigned p[4];
    enum fuzzy_state state = UNCHANGED;
    int state_set = 0;
    *cache_was_set = dir_cache_try_get_diff(a, b, &state);

    if (*cache_was_set) return state;
    else {
        /* */
        
    if (! state_set) {
        *cacheable = find_cacheable_permutation(a, alpha_packed_limits,
                                                b, alpha_packed_limits, p);
        assert(*cacheable);
        update_points_gen_params(bpar->dist[0], a, p);
        update_points_gen_params(bpar->dist[1], b, p);
        struct binomial_est_state est = get_est_state(bpar);
        state = est.state;
    }
    return state;
}
