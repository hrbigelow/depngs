/* provides memoization of classifying two dirichlets as same or
   different, based on their two largest components, and various
   confidence thresholds. */

#include "dirichlet_diff_cache.h"
#include "virtual_bound.h"
#include "yepLibrary.h"

#include <float.h>
#include <assert.h>

#define MAX(a,b) ((a) < (b) ? (b) : (a))
#define MIN(a,b) ((a) < (b) ? (a) : (b))

static struct binomial_est_params g_be_par;

unsigned alpha_packed_limits[] = {
    1<<24, 1<<20, 1<<12, 1<<8
};


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
/* UNSET must be zero */
/* enum init_phase { UNSET = 0, PENDING, SET }; */

/* struct set_flag { */
/*     pthread_mutex_t mtx; */
/*     pthread_cond_t cond; */
/* }; */

KHASH_MAP_INIT_INT64(fuzzy_hash, enum fuzzy_state);
/* KHASH_MAP_INIT_INT64(points_hash, struct counted_points); */
KHASH_MAP_INIT_INT64(bounds_hash, struct binomial_est_bounds);

struct dirichlet_diff_cache_t {
    /* unsigned n_locks; */
    /* struct set_flag *locks; */
    /* khash_t(points_hash) *dir_points; */
    /* unsigned n_vals_mtx; */
    /* pthread_mutex_t *dir_points_vals_mtx; */
    /* pthread_rwlock_t dir_points_rwlock; */
    /* struct cache_counters c; */
    /* unsigned disable_points_hash; */
    /* khash_t(bounds_hash) *bounds; */
    /* pthread_rwlock_t bounds_rwlock; */
    /* unsigned locus_count; */
    /* unsigned n_threads_bounds_active; */
    /* pthread_mutex_t locus_mtx; */
    /* pthread_cond_t locus_cond; */

} cache;


/* initializes the key from the counts, sets *packable to 1 on failure */
void 
init_alpha_packed_large(unsigned *cts, 
                        union alpha_large_key *key,
                        unsigned *packable)
{
    *packable = 
        cts[0] <= alpha_packed_limits[0]
        && cts[1] <= alpha_packed_limits[1]
        && cts[2] <= alpha_packed_limits[2]
        && cts[3] <= alpha_packed_limits[3];

    if (*packable)
        key->c = (struct alpha_packed_large){ cts[0], cts[1], cts[2], cts[3] };
}

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
                          void **reader_pars,
                          struct contig_region *qbeg,
                          struct contig_region *qend,
                          unsigned n_max_reading,
                          unsigned long max_input_mem,
                          unsigned n_threads)
{
    g_dd_par = dd_par;
    g_be_par = be_par;
    /* cache.n_locks = NUM_LOCKS; /\* this may need tuning *\/ */
    /* cache.locks = malloc(cache.n_locks * sizeof(struct set_flag)); */
    /* cache.dir_points = kh_init(points_hash); */
    /* resize so h->upper_bound is greater than max_dir_cache_items */
    /* unsigned long ub = max_dir_cache_items * (1.0 / __ac_HASH_UPPER); */
    /* kh_resize(points_hash, cache.dir_points, kroundup32(ub)); */
    /* cache.n_vals_mtx = n_threads * NUM_HASH_MTX_PER_THREAD; */
    /* cache.dir_points_vals_mtx = malloc(cache.n_vals_mtx * sizeof(pthread_mutex_t)); */
    /* pthread_rwlock_init(&cache.dir_points_rwlock, NULL); */
    /* cache.c = (struct cache_counters){  */
    /*     .n_items = 0,  */
    /*     .n_permanent_items = 0,  */
    /*     .max_items = max_dir_cache_items,  */
    /*     .n_times_cleared = 0,  */
    /*     .n_threads_points_active = n_threads,  */
    /*     .min_n_hit_to_keep = 2 */
    /* }; */

    /* pthread_cond_init(&cache.c.cond, NULL); */
    /* pthread_mutex_init(&cache.c.mtx, NULL);  */
    /* cache.b = (struct cache_counters){ 0, max_bounds_items, ULONG_MAX, 0, n_threads }; */
    /* pthread_cond_init(&cache.b.cond, NULL); */
    /* pthread_mutex_init(&cache.b.mtx, NULL); */
    /* cache.bounds = kh_init(bounds_hash); */
    /* unsigned long ub_bounds_hash =  */
    /*     (MAX_BOUNDS_PREPOP + MAX_ALPHA2_PSEUDO_PREPOP  */
    /*      + N_LOCI_TO_FREEZE) * (1.0 / __ac_HASH_UPPER); */

    /* kh_resize(bounds_hash, cache.bounds, kroundup32(ub_bounds_hash)); */

    /* pthread_rwlock_init(&cache.bounds_rwlock, NULL); */
    /* cache.locus_count = 0; */
    /* cache.n_threads_bounds_active = n_threads; */
    /* pthread_mutex_init(&cache.locus_mtx, NULL); */
    /* pthread_cond_init(&cache.locus_cond, NULL); */

    /* unsigned i; */
    /* for (i = 0; i != cache.n_locks; ++i) { */
    /*     pthread_mutex_init(&cache.locks[i].mtx, NULL); */
    /*     pthread_cond_init(&cache.locks[i].cond, NULL); */
    /* } */
    /* for (i = 0; i != cache.n_vals_mtx; ++i) */
    /*     pthread_mutex_init(&cache.dir_points_vals_mtx[i], NULL); */

    enum YepStatus status = yepLibrary_Init();
    assert(status == YepStatusOk);

    dirichlet_points_gen_init(g_dd_par.prior_alpha);

    printf("Precomputing confidence interval statistics...");
    binomial_est_init(be_par, be_par.max_sample_points, n_threads);
    printf("done.\n");

    printf("Collecting input statistics...");
    run_survey(dc_par, reader_pars, qbeg, qend, 
               n_threads, n_max_reading, max_input_mem);

    printf("done.\n");
    printf("Generating dirichlet point sets...");
    generate_point_sets(n_threads, be_par.max_sample_points);
    printf("done.\n");

    printf("Calculating bounds...");
    generate_est_bounds(n_threads);
    printf("done.\n");
}


void
dirichlet_diff_cache_free()
{
    /* unsigned i; */
    /* for (i = 0; i != cache.n_locks; ++i) { */
    /*     pthread_mutex_destroy(&cache.locks[i].mtx); */
    /*     pthread_cond_destroy(&cache.locks[i].cond); */
    /* } */
    /* free(cache.locks); */
    /* for (i = 0; i != cache.n_vals_mtx; ++i) */
    /*     pthread_mutex_destroy(&cache.dir_points_vals_mtx[i]); */

    /* free(cache.dir_points_vals_mtx); */

    /* kh_destroy(points_hash, cache.dir_points); */
    /* pthread_mutex_destroy(&cache.c.mtx); */
    /* pthread_rwlock_destroy(&cache.dir_points_rwlock); */
    /* kh_destroy(bounds_hash, cache.bounds); */
    /* pthread_rwlock_destroy(&cache.bounds_rwlock); */

    /* pthread_mutex_destroy(&cache.locus_mtx); */
    /* pthread_cond_destroy(&cache.locus_cond); */

    binomial_est_free();
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
        



/* Interpolate the interval in the beb row */
enum fuzzy_state locate_cell(struct binomial_est_bounds *beb, unsigned a1)
{
    if (beb->unchanged[0] <= a1 && a1 < beb->unchanged[1])
        return UNCHANGED;
    else if (beb->ambiguous[0] <= a1 && a1 < beb->ambiguous[1])
        return AMBIGUOUS;
    else return CHANGED;
}



/* Given a[4], b[4], and lim[4], find a permutation of {(a[p_1],
   b[p_1]), (a[p_2], b[p_2]), (a[p_3], b[p_3]), (a[p_4], b[p_4]) }
   such that both a[p_i] and b[p_i] are less than lim[p_i], if such
   permutation exists. 
   
   Modification: Given lim[2] == 1 and lim[3] == 1, The permutation of
   the last two elements does not matter.  So, we only need to specify
   p_1 and p_2 in this case.

   Also, if there is no permutation that suffices, set a flag
   appropriately.

   Assumes that lim is descending.  lim[i] >= lim[i+1]
*/
void
find_cacheable_permutation(const unsigned *a, 
                           const unsigned *b, 
                           const unsigned *lim,
                           unsigned *permutation, 
                           unsigned *perm_found)
{
    *perm_found = 1;
    /* mpi[i] (max permutation index) is the maximum position in the
       permutation (described above) that the (a, b) pair may attain
       and still stay below the limits. */
    /* min value, min_index */
    int i, j;
    int min_mpi, mpi[] = { -1, -1, -1, -1 };
    
    unsigned tot[] = { a[0] + b[0], a[1] + b[1], a[2] + b[2], a[3] + b[3] };
    
    for (i = 0; i != 4; ++i)
        for (j = 0; j != 4; ++j) {
            if (a[i] >= lim[j] || b[i] >= lim[j]) break;
            mpi[i] = j;
        }
    
    /* p[i], (the permutation of mpi), s.t. mpi[p[i]] <= mpi[p[i+1]] */

    /* mpi[j] = -1 means the value j has been placed in the
       permutation */
    for (i = 0; i != 4; ++i) {
        unsigned max_v = 0, max_i = 4;
        min_mpi = 5;
        for (j = 0; j != 4; ++j)
            if (mpi[j] >= i && mpi[j] <= min_mpi && tot[j] >= max_v) {
                max_v = tot[j];
                max_i = j;
                min_mpi = mpi[j];
            }
        if (max_i != 4)
            permutation[i] = max_i, mpi[max_i] = -1;
        else {
            *perm_found = 0;
            break;
        }
    }
}


/* estimate distance between two initialized 4D dirichlet
   distributions. Assume bpar->dist[0] and bpar->dist[1] have already
   had alpha_counts initialized, and have their size field set
   appropriately (to zero if the points buffer is invalid). Use
   caching of individual dirichlet distributions whenever possible. */
struct binomial_est_state 
get_est_state(struct bound_search_params *bpar)
{
    union alpha_large_key key[2];
    unsigned packable[2];
    unsigned *cts[2];
    unsigned i;
    khiter_t itr;
    struct points_gen_par *pgp;
    for (i = 0; i != 2; ++i) {
        pgp = bpar->dist[i]->pgen.points_gen_par;
        cts[i] = pgp->alpha_counts;
        init_alpha_packed_large(cts[i], &key[i], &packable[i]);
        struct points_buf *pb = &bpar->dist[i]->points;
        if (pb->size == 0 && packable[i]) {
            /* attempt to retrieve points from the global point sets hash */
            if ((itr = kh_get(points_h, g_points_hash, key[i].key)) 
                != kh_end(g_points_hash)) {
                pb->p = kh_val(g_points_hash, itr);
                pb->size = g_be_par.max_sample_points;
            } else {
                pb->p = pb->buf;
                pb->size = 0;
            }
        }
            /* sync_points(key[i], &bpar->dist[i]->points,  */
            /*             bpar->points_hash_frozen); */
    }

    struct binomial_est_state rval =
        binomial_quantile_est(bpar->dist[0]->pgen,
                              &bpar->dist[0]->points, 
                              bpar->dist[1]->pgen,
                              &bpar->dist[1]->points);

    /* for (i = 0; i != 2; ++i) { */
    /*     if (packable[i]) */
    /*     sync_points(key[i], &bpar->dist[i]->points,  */
    /*                 bpar->points_hash_frozen); */
    /* } */

    return rval;
}


/* The main distance function.  This will modify dist[0]'s alpha
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
    a->x = xmin, a->y = y.beta_qval_lo, a->left = NULL, a->right = b, a->down = (void *)0xDEADBEEF;

    y = pair_dist_aux(xmax, bpar);
    b->x = xmax, b->y = y.beta_qval_lo, b->left = a, b->right = NULL, b->down = (void *)0xDEADBEEF;

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
   0, 0 } and { b1, b2, 0, 0 }. */
int
get_fuzzy_state(struct bound_search_params *bpar,
                unsigned a1, unsigned a2, 
                unsigned b1, unsigned b2,
                enum fuzzy_state *state,
                unsigned *was_set)
{
    *was_set = 1;

    union bounds_key bk = { .f = { a2, b1, b2 } };
    khiter_t itr;
    struct binomial_est_bounds beb;
    int state_unset;

    if ((itr = kh_get(bounds_h, g_bounds_hash, bk.key)) != kh_end(g_bounds_hash)) {
        beb = kh_val(g_bounds_hash, itr);
        *state = locate_cell(&beb, a1);
        state_unset = 0;
    } else
        state_unset = 1;
    
    return state_unset;
}


/* test two dirichlets based on their counts. if the pattern of counts
   is cacheable, 'cacheable' is set to 1, and the cache is queried for
   the entry, and cache_was_set is set to 1 if found. */
enum fuzzy_state
cached_dirichlet_diff(unsigned *a,
                      unsigned *b,
                      struct bound_search_params *bpar,
                      unsigned *cacheable,
                      unsigned *cache_was_set)
{
    unsigned lim[] = { g_dd_par.pseudo_depth + 1, g_dd_par.pseudo_depth + 1, 1, 1 };
    unsigned p[4];
    enum fuzzy_state state = UNCHANGED;
    int state_unset;

    find_cacheable_permutation(a, b, lim, p, cacheable);
    if (*cacheable) {
        unsigned a1 = a[p[0]], a2 = a[p[1]], b1 = b[p[0]], b2 = b[p[1]];
        state_unset = get_fuzzy_state(bpar, a1, a2, b1, b2, &state, cache_was_set);
    }
    else state_unset = 1;

    if (state_unset) {
        /* not worth caching */
        unsigned lim2[] = { UINT_MAX, UINT_MAX, UINT_MAX, UINT_MAX };
        find_cacheable_permutation(a, b, lim2, p, cacheable);
        assert(*cacheable);
        update_points_gen_params(bpar->dist[0], a, p);
        update_points_gen_params(bpar->dist[1], b, p);
        struct binomial_est_state est = get_est_state(bpar);
        state = est.state;
    }
    return state;
}
