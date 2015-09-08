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
    unsigned db = dd_par.mode_batch_size, bb = be_par.batch_size;
    /* round up db to nearest multiple of bb. */
    dd_par.mode_batch_size += (db % bb) ? (bb - (db % bb)) : 0;
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

    fprintf(stdout, "%s: Caching dirichlet point sets.\n", timer_progress());
    generate_point_sets(n_threads);
    fprintf(stdout, "%s: Finished caching dirichlet point sets.\n", timer_progress());

    fprintf(stdout, "%s: Caching sample-to-REF changes.\n", timer_progress());
    generate_ref_change(n_threads);
    fprintf(stdout, "%s: Finished caching sample-to-REF changes.\n", timer_progress());

    fprintf(stdout, "%s: Caching sample-to-sample changes.\n", timer_progress());
    generate_sam_change(n_threads);
    fprintf(stdout, "%s: Finished caching sample-to-sample changes.\n", timer_progress());
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
    struct ipoint *l, *r, *d, *u;
    struct binomial_est_state est;
};

#define UNSET (void *)0xdeadbeef

/* Assume L and R point to each other horizontally. M must be a simple
   identifier.  L and R may be arbitrary expressions. */
#define INSERT_NODE_HORZ(L, M, R)               \
    do {                                        \
        struct ipoint                           \
            *lp = (L),                          \
            *rp = (R);                          \
        (M)->l = lp;                            \
        (M)->r = rp;                            \
        lp->r = (M);                            \
        rp->l = (M);                            \
    } while (0)


/* remove nd from the vertical doubly-linked list. return the new
   head. */
void
remove_vert(struct ipoint *nd, struct ipoint **hdp)
{
    struct ipoint *u = nd->u, *d = nd->d;
    nd->u = nd->d = UNSET;
    if (u) u->d = d;
    if (d) d->u = u;
    if (*hdp == nd) *hdp = d;
}


/* insert n between u and d.  u, n, and d are all pointers, n is
   valid, but u and d may not be.   */
#define INSERT_VERT_BETWEEN(d, n, u)            \
    do {                                        \
        (n)->u = (u);                           \
        (n)->d = (d);                           \
        if (u) (u)->d = (n);                    \
        if (d) (d)->u = (n);                    \
    } while (0)
        
/* insert based on the mean */
void
insert_vert_mean(struct ipoint *n, struct ipoint **hdp)
{
    if (! *hdp) {
        /* the empty list.  n becomes the head. */
        n->u = n->d = NULL;
        *hdp = n;
        return;
    }
    /* n has at least one horizontal neighbor. */
    assert(n->l || n->r);
    
    /* find a pair of vertically adjacent nodes, with at least one of
       them a horizontal neighbor of n. */
    struct ipoint *t = n->l, *u, *d;
    while (t && t->u == UNSET) t = t->l;
    if (! t) {
        t = n->r;
        while (t && t->u == UNSET) t = t->r;
    }
    assert(t);
    u = (t->u ? t->u : t);
    d = u->d;
    
    /* scan up until u is above n or becomes NULL */
    while (u && u->est.beta_mean < n->est.beta_mean)
        d = u, u = u->u;
    
    /* scan down until d is below n or becomes NULL */
    while (d && d->est.beta_mean > n->est.beta_mean)
        u = d, d = d->d;
    
    INSERT_VERT_BETWEEN(d, n, u);
    if (! u) *hdp = n;
}




/* try to insert the node into the vertical doubly-linked list, and
   update hd if necessary.  if nd cannot be inserted due to its not
   being well-ordered, return 1.  otherwise, return 0.  if n cannot be
   inserted, it means that there is one node, or two vertically
   adjacent nodes that vertically overlap n.  set n->u and n->d to
   these, or NULL if one doesn't exist.  these can then be refined and
   re-inserted. */
int
try_insert_vert(struct ipoint *n,
                struct ipoint **hdp,
                struct ipoint **d_olap,
                struct ipoint **u_olap)
{
    assert(n->u == UNSET);
    assert(n->d == UNSET);
    
    if (! *hdp) {
        /* the empty list.  n becomes the head. */
        n->u = n->d = NULL;
        *hdp = n;
        return 0;
    }
    /* n has at least one horizontal neighbor.  (n must be
       horizontally inserted before calling insert_vert) */
    assert(n->l || n->r);
    
    /* look to the left or right of n for a node that is in the
       vertical ordering. */
    struct ipoint *t = n->l, *u, *d;
    while (t && t->u == UNSET) t = t->l;
    assert(t == NULL || t->u != UNSET);
    if (! t) {
        t = n->r;
        while (t && t->u == UNSET) t = t->r;
        assert(t->u != UNSET);
    }
    assert(t);
    u = (t->u ? t->u : t);
    d = u->d;

    /* scan up until u is above n or becomes NULL */
    while (u && u->est.beta_lo < n->est.beta_hi)
        d = u, u = u->u;

    /* scan down until d is below n or becomes NULL */
    while (d && d->est.beta_hi > n->est.beta_lo)
        u = d, d = d->d;

    /* since we moved up and then down, we re-test u to see if it is
       still a non-overlapping upper neighbor. */
    if (u == NULL || u->est.beta_lo > n->est.beta_hi) {
        INSERT_VERT_BETWEEN(d, n, u);
        if (! u) *hdp = n;
        return 0;
    } else {
        /* n hasn't been inserted vertically, but set n->u and n->d to
           the nodes that are the best candidates */
        *u_olap = (u && u->est.beta_lo < n->est.beta_hi) ? u : NULL;
        *d_olap = (d && n->est.beta_lo < d->est.beta_hi) ? d : NULL;
        return 1;
    }
}


/* recursively try to insert n vertically.  if it fails, produce
   additional samples to narrow the error bars, and repeat. augments
   with new samples not only the node being inserted but the one or
   two nodes that overlap with it vertically.  insert using the mean
   if resampling is exhausted and there is still overlap. */
void
insert_vert(struct ipoint *q, struct ipoint **hdp,
            const unsigned *bounds, unsigned n_add,
            unsigned n_max_trials) {
    
    unsigned perm_alpha1[] = { 0, bounds[0], 0, 0 };
    unsigned perm_alpha2[] = { bounds[1], bounds[2], 0, 0 };
    POINT *buf = malloc(sizeof(POINT) * n_add * 2);
    POINT *buf1 = buf, *buf2 = buf + n_add;
    struct ipoint *d, *u;
    while (try_insert_vert(q, hdp, &d, &u)) {
        if (u && u->est.n_trials < n_max_trials) {
            remove_vert(u, hdp);
            perm_alpha1[0] = u->x;
            gen_dir_points(perm_alpha1, buf1, n_add);
            gen_dir_points(perm_alpha2, buf2, n_add);
            add_binomial_trials(buf1, buf2, n_add, &u->est);
            insert_vert(u, hdp, bounds, n_add, n_max_trials);
        } else if (d && d->est.n_trials < n_max_trials) {
            remove_vert(d, hdp);
            perm_alpha1[0] = d->x;
            gen_dir_points(perm_alpha1, buf1, n_add);
            gen_dir_points(perm_alpha2, buf2, n_add);
            add_binomial_trials(buf1, buf2, n_add, &d->est);
            insert_vert(d, hdp, bounds, n_add, n_max_trials);
        } else if (q->est.n_trials < n_max_trials) {
            perm_alpha1[0] = q->x;
            gen_dir_points(perm_alpha1, buf1, n_add);
            gen_dir_points(perm_alpha2, buf2, n_add);
            add_binomial_trials(buf1, buf2, n_add, &q->est);
        }
        else {
            insert_vert_mean(q, hdp);
            break;
        }
    }
    free(buf);
}


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
        ? (est.beta_lo < bsp->query_beta ? 1 : 0)
        : (est.beta_hi < bsp->query_beta ? 1 : 0);
}


/* Number of highest points to explore by flanking, used in
   noisy_mode. */
#define MAX_NEW_NODES 1000

/* find the mode in a semi-robust way.  for u and d members,
   0xDEADBEEF is used to mean 'uninitialized', while NULL means 'last
   in the chain' */
static unsigned
noisy_mode(unsigned xmin, unsigned xend, unsigned n_add,
           unsigned n_max_trials, void *par)
{
    /* create the two anchors.  (we malloc these because they become
       part of the list, so it is easier to free the whole thing at
       the end) */
    struct ipoint *hd = malloc(sizeof(struct ipoint));
    struct ipoint *nd = malloc(sizeof(struct ipoint));
    struct bound_search_params *bsp = par;
    
    unsigned xmax = xend - 1;

    /* create a pair of nodes at each end, horizontally ordered. hd->u
       and hd->d are set to NULL to indicate that hd is vertically
       ordered in a single-node list. */
    *hd = (struct ipoint){
        .est = pair_dist_aux(bsp->dist[0], bsp->dist[1], xmin),
        .x = xmin, .l = NULL, .r = nd,
        .u = NULL, .d = NULL
    };
    
    *nd = (struct ipoint){
        .est = pair_dist_aux(bsp->dist[0], bsp->dist[1], xmax),
        .x = xmax, .l = hd, .r = NULL,
        .u = UNSET, .d = UNSET
    };
    
    unsigned bounds[] = { bsp->dist[0]->perm_alpha[1],
                          bsp->dist[1]->perm_alpha[0],
                          bsp->dist[1]->perm_alpha[1] };
    
    insert_vert(nd, &hd, bounds, n_add, n_max_trials);
    
    struct ipoint *newnodes[MAX_NEW_NODES];
    
    while (1) {
        /* create flanking nodes for all nodes that overlap each other
           vertically. stop when the next lower node does not
           overlap. */
        struct ipoint *cen = hd;
        unsigned i = 0, x;
        while (i < MAX_NEW_NODES - 1) {
            /* find a position halfway between cen and cen->l,
               horizontally insert a new node there if space
               exists. */
            if (cen->l
                && (x = (cen->x + cen->l->x) / 2) != cen->x
                && x != cen->l->x) {
                struct ipoint *nn = malloc(sizeof(struct ipoint));
                nn->x = x;
                nn->est = pair_dist_aux(bsp->dist[0], bsp->dist[1], nn->x);
                INSERT_NODE_HORZ(cen->l, nn, cen);
                nn->u = nn->d = UNSET;
                newnodes[i++] = nn;
            }
            /* find a position halfway between cen and cen->r,
               horizontally insert a new node there if space
               permits. */
            if (cen->r
                && (x = (cen->x + cen->r->x) / 2) != cen->x
                && x != cen->r->x) {
                struct ipoint *nn = malloc(sizeof(struct ipoint));
                nn->x = x;
                nn->est = pair_dist_aux(bsp->dist[0], bsp->dist[1], nn->x);
                INSERT_NODE_HORZ(cen, nn, cen->r);
                nn->u = nn->d = UNSET;
                newnodes[i++] = nn;
            }
            /* we break when there is a non-overlapping break in the
               vertical intervals. */
            if (! cen->d || cen->est.beta_lo > cen->d->est.beta_hi)
                break;
            cen = cen->d;
        }
        /* Now, up to 4 new nodes are created.  They will be the
           immediate neighbors of the top two nodes */
        unsigned n_new = i;
        for (i = 0; i != n_new; ++i)
            insert_vert(newnodes[i], &hd, bounds, n_add, n_max_trials);
        
        /* We want these to represent distance to neighbor along x
           axis. Where there is no neighbor, consider it to be '1' */
        unsigned lx = (int)hd->x - (hd->l ? hd->l->x : -1);
        unsigned rx = (hd->r ? hd->r->x : xend) - hd->x;
        
        /* stop iff both neighbors are as close as possible */
        if (lx == 1 && rx == 1) break;
    }
    unsigned topx = hd->x;
    
    /* print out nodes from left to right */
    struct ipoint *tp = hd;
    while (tp->l) tp = tp->l;
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
                tp->est.beta_lo, tp->est.beta_hi);
        tp = tp->r;
    }

    /* clean up */
    struct ipoint *p;
    while (hd) {
        p = hd;
        hd = hd->d;
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
    unsigned xmode = noisy_mode(0,
                                g_dd_par.xmax,
                                g_dd_par.mode_batch_size,
                                g_dd_par.max_bernoulli_trials,
                                bsp);

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
