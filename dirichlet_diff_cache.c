/* provides memoization of classifying two dirichlets as same or
   different, based on their two largest components, and various
   confidence thresholds. */

#include "binomial_est.h"
#include "virtual_bound.h"
#include "dirichlet_diff_cache.h"
#include "dirichlet_points_gen.h"
#include "khash.h"
#include "cache.h"

#include <math.h>
#include <string.h>
#include <pthread.h>
#include <float.h>
#include <inttypes.h>
#include <assert.h>

#define MAX(a,b) ((a) < (b) ? (b) : (a))
#define MIN(a,b) ((a) < (b) ? (a) : (b))



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
enum init_phase { UNSET = 0, PENDING, SET };

struct binomial_est_bounds {
    enum init_phase state;
    int32_t ambiguous[2];
    int32_t unchanged[2];
};

struct set_flag {
    pthread_mutex_t mtx;
    pthread_cond_t cond;
};

struct counted_points {
    struct points_buf pts;
    unsigned hitcount;
};

struct counted_bounds {
    struct binomial_est_bounds beb;
    unsigned hitcount;
};

KHASH_MAP_INIT_INT64(fuzzy_hash, enum fuzzy_state)
KHASH_MAP_INIT_INT64(points_hash, struct counted_points);
KHASH_MAP_INIT_INT64(bounds_hash, struct counted_bounds);

struct cache_counters {
    unsigned long n_items, max_items;
    unsigned long n_last_cleared_entries;
    unsigned n_times_cleared;
    unsigned n_unset_flags;
    pthread_cond_t cond;
    pthread_mutex_t mtx;
};

struct dirichlet_diff_cache_t {
    unsigned max_depth;
    unsigned batch_size;
    double post_confidence;
    double beta_confidence;
    double min_dirichlet_dist;
    unsigned max_sample_points;
    unsigned n_locks;
    struct set_flag *locks;
    khash_t(points_hash) *dir_points;
    pthread_mutex_t dir_points_mtx;
    unsigned n_vals_mtx;
    pthread_mutex_t *dir_points_vals_mtx;
    struct cache_counters c;
    khash_t(bounds_hash) *bounds;
    pthread_mutex_t bounds_mtx;
    pthread_cond_t bounds_cond;
    struct cache_counters b;

    /* struct binomial_est_bounds *buf1, **buf2, ***bounds; */
    /* unsigned long max_secondary_cache_size; */
    /* khash_t(fuzzy_hash) *hash; */
    /* pthread_mutex_t hash_mtx; */
} cache;


struct {
    unsigned long hit, miss;
} cache_stats;


/* Arbitrary value when it makes sense to clear entries from the
   cache. */

/* The first thread to get to 'do_freeze' waits, letting subsequent
   threads enter this function and wait.  The last thread entering
   will decrement n_unset_flags down to zero, and then broadcast to
   others. */
unsigned hash_frozen_aux(struct cache_counters *cc,
                         unsigned long min_last_cleared_entries,
                         unsigned long max_times_cleared)
{
    pthread_mutex_lock(&cc->mtx);
    unsigned do_freeze = 
        cc->n_last_cleared_entries <= min_last_cleared_entries
        || cc->n_times_cleared > max_times_cleared;
    if (do_freeze) 
    {
        cc->n_unset_flags--;
        if (cc->n_unset_flags)
            pthread_cond_wait(&cc->cond, &cc->mtx);
        else
            pthread_cond_broadcast(&cc->cond);
    }
    pthread_mutex_unlock(&cc->mtx);
    return do_freeze;
}


#define MIN_LAST_CLEARED_POINTS_ENTRIES 1000
#define MAX_TIMES_CLEARED 50
unsigned points_hash_frozen()
{
    return hash_frozen_aux(&cache.c, 
                           MIN_LAST_CLEARED_POINTS_ENTRIES,
                           MAX_TIMES_CLEARED);
}

#define MIN_LAST_CLEARED_BOUNDS_ENTRIES 100
unsigned bounds_hash_frozen()
{
    return hash_frozen_aux(&cache.b,
                           MIN_LAST_CLEARED_BOUNDS_ENTRIES,
                           MAX_TIMES_CLEARED);
}

void print_cache_stats()
{
    pthread_mutex_lock(&cache.dir_points_mtx);
    time_t cal = time(NULL);
    char *ts = strdup(ctime(&cal));
    ts[strlen(ts)-1] = '\0';
    fprintf(stderr, "%s: Dirichlet cache: %lu hits, %lu misses, %5.3f hit rate\n",
            ts,
            cache_stats.hit, cache_stats.miss,
            100.0 * (float)cache_stats.hit 
            / (float)(cache_stats.hit + cache_stats.miss));
    cache_stats.hit = 0;
    cache_stats.miss = 0;
    free(ts);
    pthread_mutex_unlock(&cache.dir_points_mtx);
}


struct alpha_packed_large {
    unsigned a0 :24;
    unsigned a1 :20;
    unsigned a2 :12;
    unsigned a3 :8;
};

union alpha_large_key {
    struct alpha_packed_large c;
    uint64_t raw;
};

/* initializes the key from the counts, sets ret to 1 on failure */
void init_alpha_packed_large(unsigned *cts, union alpha_large_key *key,
                             unsigned *packable)
{
    *packable = cts[0] <= (1<<24) 
        && cts[1] <= (1<<20) 
        && cts[2] <= (1<<12) 
        && cts[3] <= (1<<8);

    if (*packable)
        key->c = (struct alpha_packed_large){ cts[0], cts[1], cts[2], cts[3] };
}

struct alpha_packed {
    unsigned a0 :12; /* 4096 */
    unsigned a1 :10; /* 1024 */
    unsigned a2 :6;  /* 64 */
    unsigned a3 :4;  /* 16 */
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
        

#define NUM_LOCKS 10
#define NUM_HASH_MTX_PER_THREAD 1

/* bounds_cache[a2][b2][b1].  sets each element to have 'state' = UNSET */
void dirichlet_diff_init(unsigned max_depth,
                         unsigned batch_size,
                         double post_confidence,
                         double beta_confidence,
                         double min_dirichlet_dist,
                         unsigned max_sample_points,
                         unsigned long max_dir_cache_items,
                         unsigned long max_bounds_items,
                         unsigned n_threads)
{
    cache.max_depth = max_depth;
    cache.batch_size = batch_size;
    cache.post_confidence = post_confidence;
    cache.beta_confidence = beta_confidence;
    cache.min_dirichlet_dist = min_dirichlet_dist;
    cache.max_sample_points = max_sample_points;
    cache.n_locks = NUM_LOCKS; /* this may need tuning */
    cache.locks = malloc(cache.n_locks * sizeof(struct set_flag));
    cache.dir_points = kh_init(points_hash);
    /* resize so h->upper_bound is greater than max_dir_cache_items */
    unsigned long ub = max_dir_cache_items * (1.0 / __ac_HASH_UPPER);
    kh_resize(points_hash, cache.dir_points, kroundup32(ub));
    pthread_mutex_init(&cache.dir_points_mtx, NULL);
    cache.n_vals_mtx = n_threads * NUM_HASH_MTX_PER_THREAD;
    cache.dir_points_vals_mtx = malloc(cache.n_vals_mtx * sizeof(pthread_mutex_t));

    cache.c = (struct cache_counters){ 0, max_dir_cache_items, ULONG_MAX, 0, n_threads };
    pthread_cond_init(&cache.c.cond, NULL);
    pthread_mutex_init(&cache.c.mtx, NULL); 
    cache.b = (struct cache_counters){ 0, max_bounds_items, ULONG_MAX, 0, n_threads };
    pthread_cond_init(&cache.b.cond, NULL);
    pthread_mutex_init(&cache.b.mtx, NULL);
    cache.bounds = kh_init(bounds_hash);
    pthread_mutex_init(&cache.bounds_mtx, NULL);
    pthread_cond_init(&cache.bounds_cond, NULL);

    /* cache.buf1 = calloc(max2 * max2 * max1, sizeof(struct binomial_est_bounds)); */
    /* cache.buf2 = malloc(max2 * max2 * sizeof(struct binomial_est_bounds *)); */
    /* cache.bounds = malloc(max2 * sizeof(struct binomial_est_bounds **)); */
    /* cache.max_secondary_cache_size = max_secondary_cache_size; */
    /* cache.hash = kh_init(fuzzy_hash); */
    /* pthread_mutex_init(&cache.hash_mtx, NULL); */

    /* unsigned i, j, k, l; */
    /* unsigned M = max2 * max2; */
    /* for (i = 0, j = 0; i != M; ++i, j += max1) */
    /*     cache.buf2[i] = cache.buf1 + j; */

    /* for (k = 0, l = 0; k != max2; ++k, l += max2) */
    /*     cache.bounds[k] = cache.buf2 + l; */

    unsigned i;
    for (i = 0; i != cache.n_locks; ++i)
    {
        pthread_mutex_init(&cache.locks[i].mtx, NULL);
        pthread_cond_init(&cache.locks[i].cond, NULL);
    }
    for (i = 0; i != cache.n_vals_mtx; ++i)
        pthread_mutex_init(&cache.dir_points_vals_mtx[i], NULL);

    cache_stats.hit = 0;
    cache_stats.miss = 0;
}


void dirichlet_diff_free()
{
    /* free(cache.bounds); */
    /* free(cache.buf2); */
    /* free(cache.buf1); */

    unsigned i;
    for (i = 0; i != cache.n_locks; ++i)
    {
        pthread_mutex_destroy(&cache.locks[i].mtx);
        pthread_cond_destroy(&cache.locks[i].cond);
    }
    free(cache.locks);
    for (i = 0; i != cache.n_vals_mtx; ++i)
        pthread_mutex_destroy(&cache.dir_points_vals_mtx[i]);

    free(cache.dir_points_vals_mtx);

    kh_destroy(points_hash, cache.dir_points);
    pthread_mutex_destroy(&cache.c.mtx);
    pthread_mutex_destroy(&cache.dir_points_mtx);
    kh_destroy(bounds_hash, cache.bounds);
    pthread_mutex_destroy(&cache.bounds_mtx);
    pthread_cond_destroy(&cache.bounds_cond);

    /* kh_destroy(fuzzy_hash, cache.hash); */
    /* pthread_mutex_destroy(&cache.hash_mtx); */
}


void alloc_distrib_points(struct distrib_points *dpts)
{
    unsigned msp = cache.max_sample_points;
    dpts->pgen = (struct points_gen){ 
        malloc(sizeof(struct points_gen_par)),
        gen_dirichlet_points_wrapper, 
        calc_post_to_dir_ratio
    };
    ((struct points_gen_par *)dpts->pgen.points_gen_par)->randgen = 
        gsl_rng_alloc(gsl_rng_taus);
    dpts->points = (struct points_buf){ (POINT *)malloc(sizeof(POINT) * msp), 0, msp };
    dpts->weights = (struct weights_buf){ (double *)malloc(sizeof(double) * msp), 0, msp };
}


void free_distrib_points(struct distrib_points *dpts)
{
    gsl_rng_free(((struct points_gen_par *)dpts->pgen.points_gen_par)->randgen);
    free((struct points_gen_par *)dpts->pgen.points_gen_par);
    free(dpts->points.buf);
    free(dpts->weights.buf);
}


/* update dpts parameters and discard existing points.  if alpha is
   identical to existing parameters, do nothing. */
void update_points_gen_params(struct distrib_points *dpts,
                              unsigned *alpha_counts,
                              unsigned *perm)
{
    struct points_gen_par *dp = dpts->pgen.points_gen_par;
    unsigned i, change = 0;
    for (i = 0; i != NUM_NUCS; ++i)
    {
        if (dp->alpha_counts[i] != alpha_counts[perm[i]]
            || dp->alpha_perm[i] != perm[i]) change = 1;
        dp->alpha_counts[i] = alpha_counts[perm[i]];
        dp->alpha_perm[i] = perm[i];
    }

    if (change)
    {
        dpts->points.size = 0;
        dpts->weights.size = 0;
    }
}

/* set a single alpha, and discard points if there is a change */
void set_dirichlet_alpha_single(struct distrib_points *dpts, unsigned i, unsigned v)
{
    struct points_gen_par *dp = dpts->pgen.points_gen_par;
    if (dp->alpha_counts[i] != v)
    {
        dp->alpha_counts[i] = v;
        dpts->points.size = 0;
        dpts->weights.size = 0;
    }
}

void read_diststats_line(FILE *fh, 
                         unsigned *a2,
                         unsigned *b1,
                         unsigned *b2,
                         struct binomial_est_bounds *beb)
{
    int n;
    n = fscanf(fh, 
               "%u\t%u\t%u\t%"SCNd32"\t%"SCNd32"\t%"SCNd32"\t%"SCNd32"\n", 
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
/*
void parse_diststats_header(FILE *diststats_fh, double *prior_alpha)
{
    int n;
    n = fscanf(diststats_fh, 
               "# max1: %u\n"
               "# max2: %u\n"
               "# max_sample_points: %u\n"
               "# min_dirichlet_dist: %lf\n"
               "# post_confidence: %lf\n"
               "# beta_confidence: %lf\n"
               "# prior_alpha: %lf\n",
               &cache.max1, 
               &cache.max2, 
               &cache.max_sample_points,
               &cache.min_dirichlet_dist,
               &cache.post_confidence,
               &cache.beta_confidence,
               prior_alpha);
    if (n != 7)
    {
        fprintf(stderr, "error parsing diststats header at %s: %ul\n", __FILE__, __LINE__);
        exit(1);
    }
               
}
*/

 /*
void write_diststats_header(FILE *diststats_fh)
{
    fprintf(diststats_fh, 
            "# max1: %u\n"
            "# max2: %u\n"
            "# max_sample_points: %u\n"
            "# min_dirichlet_dist: %lf\n"
            "# post_confidence: %lf\n"
            "# beta_confidence: %lf\n"
            "# prior_alpha: %g\n",
            cache.max1, 
            cache.max2, 
            cache.max_sample_points,
            cache.min_dirichlet_dist,
            cache.post_confidence,
            cache.beta_confidence,
            get_alpha_prior());
}
 */

/* initialize internal bounds_cache */
/*
void parse_diststats_body(FILE *diststats_fh, unsigned max1, unsigned max2)
{
    struct binomial_est_bounds beb;
    unsigned a2, b1, b2;
    while (! feof(diststats_fh))
    {
        read_diststats_line(diststats_fh, &a2, &b1, &b2, &beb);
        if (a2 < max2 && b2 < max2 && b1 < max1)
            cache.bounds[a2][b2][b1] = beb;
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
*/

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
void find_cacheable_permutation(const unsigned *a, 
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
        for (j = 0; j != 4; ++j)
        {
            if (a[i] >= lim[j] || b[i] >= lim[j]) break;
            mpi[i] = j;
        }
    
    /* p[i], (the permutation of mpi), s.t. mpi[p[i]] <= mpi[p[i+1]] */

    /* mpi[j] = -1 means the value j has been placed in the
       permutation */
    for (i = 0; i != 4; ++i)
    {
        unsigned max_v = 0, max_i = 4;
        min_mpi = 5;
        for (j = 0; j != 4; ++j)
            if (mpi[j] >= i && mpi[j] <= min_mpi && tot[j] >= max_v)
            {
                max_v = tot[j];
                max_i = j;
                min_mpi = mpi[j];
            }
        if (max_i != 4)
            permutation[i] = max_i, mpi[max_i] = -1;
        else
        {
            *perm_found = 0;
            break;
        }
    }
}


/* simple hash function for distributing the dir_points_vals_mtx in
   use. */
unsigned dir_points_key_group(union alpha_large_key key)
{
    uint32_t h = ((key.raw)>>33^(key.raw)^(key.raw)<<11);
    return h % cache.n_vals_mtx;
}


/* append any points in src beyond dest->size, growing dest as
   necessary. */
void append_extra_points(struct points_buf *dest, struct points_buf *src)
{
    assert(src->size >= dest->size);
    unsigned end = dest->size, add = src->size - end;
    dest->size = src->size;
    ALLOC_GROW(dest->buf, src->size, dest->alloc);
    memcpy(dest->buf + end, src->buf + end, add * sizeof(dest->buf[0]));
}

/* initialize points from cache.  we only call this if the desired set
   of points is not populated.  once points_hash_frozen is true, you
   do not need any mutex to do hash lookup operations.  however, to
   actually read the contents of the value in the hash, you
   do. finally, if val.size == cache.max_sample_points, we know that
   val.buf will remain constant as well. */
void safe_get_cached_points(union alpha_large_key key,
                            struct points_buf *points,
                            unsigned points_hash_frozen)
{
    assert(points->size == 0);
    khiter_t itr;

    pthread_mutex_t *grp_mtx = 
        &cache.dir_points_vals_mtx[dir_points_key_group(key)];

    if (points_hash_frozen)
    {
        /* no need for a mutex to query the hash */
        itr = kh_get(points_hash, cache.dir_points, key.raw);
        if (itr != kh_end(cache.dir_points))
        {
            pthread_mutex_lock(grp_mtx);
            struct counted_points stored = kh_value(cache.dir_points, itr);
            if (stored.pts.size == cache.max_sample_points)
            {
                /* the stored buffer is full and so its contents won't
                   change. */
                pthread_mutex_unlock(grp_mtx);
                append_extra_points(points, &stored.pts);
                pthread_mutex_lock(grp_mtx);
                stored.hitcount++;
                kh_value(cache.dir_points, itr) = stored;
                pthread_mutex_unlock(grp_mtx);
            }
            else
                pthread_mutex_unlock(grp_mtx);
        }
        else 
            ; /* nothing to do. */
    }
    else
    {
        pthread_mutex_lock(grp_mtx);
        itr = kh_get(points_hash, cache.dir_points, key.raw);
        if (itr != kh_end(cache.dir_points))
        {
            struct counted_points stored = kh_value(cache.dir_points, itr);
            append_extra_points(points, &stored.pts);
            stored.hitcount++;
            kh_value(cache.dir_points, itr) = stored;
            pthread_mutex_unlock(grp_mtx);
            /* ++cache_stats.hit; */
        }
        else
        {
            /* ++cache_stats.miss; */
            pthread_mutex_unlock(grp_mtx);
            /* nothing to do */
        }
    }
}


/* clear entries if necessary */
void winnow_points_hash(unsigned max_hitcount)
{
    pthread_mutex_lock(&cache.c.mtx);

    if (cache.c.n_items >= cache.c.max_items
        && cache.c.n_last_cleared_entries > MIN_LAST_CLEARED_POINTS_ENTRIES
        && cache.c.n_times_cleared < MAX_TIMES_CLEARED)
    {
        struct counted_points val;
        khiter_t itr;

        pthread_mutex_lock(&cache.dir_points_mtx);
        khint_t old_size = kh_size(cache.dir_points);
        for (itr = kh_begin(cache.dir_points);
             itr != kh_end(cache.dir_points); ++itr)
        {
            if (kh_exist(cache.dir_points, itr))
            {
                val = kh_value(cache.dir_points, itr);
                if (val.hitcount <= max_hitcount)
                {
                    kh_del(points_hash, cache.dir_points, itr);
                    if (val.pts.buf) free(val.pts.buf);
                }
            }
        }
        khint_t new_size = kh_size(cache.dir_points);
        pthread_mutex_unlock(&cache.dir_points_mtx);

        cache.c.n_items = new_size;
        cache.c.n_last_cleared_entries = old_size - new_size;
        ++cache.c.n_times_cleared;

        time_t cal = time(NULL);
        char *ts = strdup(ctime(&cal));
        ts[strlen(ts)-1] = '\0';
        fprintf(stderr, "%s: cache.n_last_cleared_entries = %lu\n",
                ts, cache.c.n_last_cleared_entries);
        free(ts);
    }

    pthread_mutex_unlock(&cache.c.mtx);
}


/* remove low hit entries from the bounds hash */
void winnow_bounds_hash(unsigned max_hitcount)
{
    pthread_mutex_lock(&cache.b.mtx);

    if (cache.b.n_items >= cache.b.max_items
        && cache.b.n_last_cleared_entries > MIN_LAST_CLEARED_BOUNDS_ENTRIES
        && cache.b.n_times_cleared < MAX_TIMES_CLEARED)
    {
        struct counted_bounds cb;
        khiter_t itr;

        pthread_mutex_lock(&cache.bounds_mtx);
        khint_t old_size = kh_size(cache.bounds);
        for (itr = kh_begin(cache.bounds);
             itr != kh_end(cache.bounds); ++itr)
        {
            if (kh_exist(cache.bounds, itr))
            {
                cb = kh_value(cache.bounds, itr);
                if (cb.hitcount <= max_hitcount)
                    kh_del(bounds_hash, cache.bounds, itr);
            }
        }
        khint_t new_size = kh_size(cache.bounds);
        pthread_mutex_unlock(&cache.bounds_mtx);

        cache.b.n_items = new_size;
        cache.b.n_last_cleared_entries = old_size - new_size;
        ++cache.b.n_times_cleared;

        time_t cal = time(NULL);
        char *ts = strdup(ctime(&cal));
        ts[strlen(ts)-1] = '\0';
        fprintf(stderr, "%s: cache.b.n_last_cleared_entries = %lu\n",
                ts, cache.b.n_last_cleared_entries);
        free(ts);
    }

    pthread_mutex_unlock(&cache.b.mtx);
}


/* store points in cache.  if points->size > stored_size, replace hash
   value with a copy of points.  outside of the mutex, hash values are
   guaranteed to be constant. a missing value is deemed to have a
   stored_size of 0. 

   the mutex logic is as follows: at a top level, if
   points_hash_frozen = 1, we know that cache.dir_points will not have
   any more keys added or deleted.  only kh_get, kh_end, and kh_val
   are called in this situation.   */
void safe_store_cached_points(union alpha_large_key key,
                              struct points_buf *points,
                              unsigned points_hash_frozen)
{
    assert(points->size > 0);
    khiter_t itr;
    pthread_mutex_t *grp_mtx = 
        &cache.dir_points_vals_mtx[dir_points_key_group(key)];

    if (points_hash_frozen)
    {
        itr = kh_get(points_hash, cache.dir_points, key.raw);
        if (itr != kh_end(cache.dir_points))
        {
            pthread_mutex_lock(grp_mtx);
            struct counted_points stored = kh_val(cache.dir_points, itr);
            if (stored.pts.size < points->size)
            {
                /* copy additional points over to stored */
                append_extra_points(&stored.pts, points);
                kh_val(cache.dir_points, itr) = stored;
                pthread_mutex_unlock(grp_mtx);
                
            }
            else
                pthread_mutex_unlock(grp_mtx);
        }
        else
            ; /* do nothing */
    }
    else
    {
        /* hash is not frozen.  try to find the key */
        struct counted_points stored;
        pthread_mutex_lock(grp_mtx);
        itr = kh_get(points_hash, cache.dir_points, key.raw);
        if (itr != kh_end(cache.dir_points))
        {
            stored = kh_val(cache.dir_points, itr);
            pthread_mutex_unlock(grp_mtx);

            if (stored.pts.size < points->size)
            {
                append_extra_points(&stored.pts, points);
                kh_val(cache.dir_points, itr) = stored;
            }                
            pthread_mutex_unlock(grp_mtx);
        }
        else
        {
            /* key did not exist */
            pthread_mutex_unlock(grp_mtx);
            struct counted_points new_buf = {
                .pts = { malloc(points->alloc * sizeof(points->buf[0])), points->size, points->alloc },
                .hitcount = 1
            };
            memcpy(new_buf.pts.buf, points->buf, points->size * sizeof(points->buf[0]));

            pthread_mutex_lock(&cache.dir_points_mtx);
            unsigned i = 0, n_locks_needed = cache.n_vals_mtx;
            while (n_locks_needed)
            {
                int rv = pthread_mutex_trylock(&cache.dir_points_vals_mtx[i]);
                if (rv == 0) --n_locks_needed;
                ++i; i %= cache.n_vals_mtx;
            }

            /* test again while safely locked */
            itr = kh_get(points_hash, cache.dir_points, key.raw);
            if (itr == kh_end(cache.dir_points))
            {            
                int ret;
                itr = kh_put(points_hash, cache.dir_points, key.raw, &ret);
                assert(ret == 1 || ret == 2); /* key did not previously exist */
                kh_val(cache.dir_points, itr) = new_buf;
            }
            else
                free(new_buf.pts.buf);

            for (i = 0; i != cache.n_vals_mtx; ++i)
                (void)pthread_mutex_unlock(&cache.dir_points_vals_mtx[i]);
            pthread_mutex_unlock(&cache.dir_points_mtx);
            
            pthread_mutex_lock(&cache.c.mtx);
            cache.c.n_items = kh_size(cache.dir_points);
            pthread_mutex_unlock(&cache.c.mtx);
        }
    }

    /* clean out rarely used entries */
    if (! points_hash_frozen) winnow_points_hash(1);
}


/* calculate the distance between two initialized
   distributions. Assume bpar->dist[0] and bpar->dist[1] have already
   had alpha_counts initialized, and have their size field set
   appropriately (to zero if the points buffer is invalid). Use
   caching of individual dirichlet distributions whenever possible. */
struct binomial_est_state 
get_est_state(struct binomial_est_params *bpar)
{
    union alpha_large_key a_key, b_key;
    unsigned a_packable, b_packable;
    unsigned *a_cts = ((struct points_gen_par *)bpar->dist[0]->pgen.points_gen_par)->alpha_counts;
    unsigned *b_cts = ((struct points_gen_par *)bpar->dist[1]->pgen.points_gen_par)->alpha_counts;
    init_alpha_packed_large(a_cts, &a_key, &a_packable);
    init_alpha_packed_large(b_cts, &b_key, &b_packable);

    if (bpar->dist[0]->points.size == 0 && a_packable)
        safe_get_cached_points(a_key, &bpar->dist[0]->points, bpar->points_hash_frozen);

    if (bpar->dist[1]->points.size == 0 && b_packable)
        safe_get_cached_points(b_key, &bpar->dist[1]->points, bpar->points_hash_frozen);

    struct binomial_est_state rval =
        binomial_quantile_est(cache.max_sample_points, 
                              cache.min_dirichlet_dist, 
                              cache.post_confidence,
                              cache.beta_confidence, 
                              bpar->dist[0]->pgen,
                              &bpar->dist[0]->points, 
                              bpar->dist[1]->pgen,
                              &bpar->dist[1]->points,
                              cache.batch_size);

    if (a_packable) 
        safe_store_cached_points(a_key, &bpar->dist[0]->points, bpar->points_hash_frozen);
    if (b_packable)
        safe_store_cached_points(b_key, &bpar->dist[1]->points, bpar->points_hash_frozen);

    return rval;
}


/* The main distance function.  This will modify dist[0]'s alpha
   parameters and thus discard its points. It will assume all 7 other
   alpha components are set as intended.  */
struct binomial_est_state pair_dist_aux(unsigned a1, void *par)
{
    struct binomial_est_params *bpar = par;
    set_dirichlet_alpha_single(bpar->dist[0], 0, a1);
    return get_est_state(bpar);
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
    /* struct points_gen_par  */
    /*     *d0 = bb->dist[0]->pgen.points_gen_par, */
    /*     *d1 = bb->dist[1]->pgen.points_gen_par; */

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



/* For two dirichlet distributions A = { x+p, a2+p, p, p } and B = {
   b1+p, b2+p, p, p }, virtually find the values of x in [0, max1)
   that denote the intervals where A and B are UNCHANGED, and where
   they are AMBIGUOUS.  The expected underlying pattern as a function
   of x is some number of C, AC, A, AU, U, AU, A, AC, C.  (see enum
   fuzzy_state in binomial_est.h) */
void initialize_est_bounds(unsigned a2, unsigned b1, unsigned b2,
                           struct binomial_est_params *bpar,
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
    unsigned xmode = noisy_mode(0, cache.max_depth, bpar);

    bpar->use_low_beta = 0;
    bpar->query_beta = 1.0 - cache.post_confidence;
    beb->ambiguous[0] = virtual_lower_bound(0, xmode, elem_is_less, bpar);
    beb->ambiguous[1] = virtual_upper_bound(xmode, cache.max_depth, elem_is_less, bpar);

    bpar->use_low_beta = 1;
    bpar->query_beta = cache.post_confidence;
    beb->unchanged[0] = virtual_lower_bound(beb->ambiguous[0], xmode, elem_is_less, bpar);
    beb->unchanged[1] = virtual_upper_bound(xmode, beb->ambiguous[1], elem_is_less, bpar);
    
}


union bounds_key {
    struct {
        unsigned a2:20;
        unsigned b1:32;
        unsigned b2:20;
    } f;
    int64_t val;
};

/* find and optionally cache a distance classification of { a1, a2, 0,
   0 } and { b1, b2, 0, 0 } (or an equivalent tandem
   permutation). returns 0 on success, or 1 if no state could be
   set. */
int
get_fuzzy_state(struct binomial_est_params *bpar,
                unsigned a1, unsigned a2, 
                unsigned b1, unsigned b2,
                enum fuzzy_state *state,
                unsigned hash_frozen,
                unsigned *was_set)
{
    *was_set = 1;

    union bounds_key bk = { .f = { a2, b1, b2 } };
    khiter_t k;
    struct counted_bounds cb;
    int state_unset;

    if (hash_frozen)
    {
        k = kh_get(bounds_hash, cache.bounds, bk.val);
        if (k != kh_end(cache.bounds))
        {
            cb = kh_val(cache.bounds, k);
            *state = locate_cell(&cb.beb, a1);
            state_unset = 0;
        }
        else
        {
            *was_set = 0;
            state_unset = 1;
        }
    }
    else
    {
        pthread_mutex_lock(&cache.bounds_mtx);
        k = kh_get(bounds_hash, cache.bounds, bk.val);
        if (k == kh_end(cache.bounds))
        {
            /* key not found.  */
            *was_set = 0;
            int ret;
            k = kh_put(bounds_hash, cache.bounds, bk.val, &ret);
            cb.hitcount = 1;
            cb.beb.state = PENDING;
            kh_val(cache.bounds, k) = cb;

            pthread_mutex_lock(&cache.b.mtx);
            cache.b.n_items = kh_size(cache.bounds);
            pthread_mutex_unlock(&cache.b.mtx);

            pthread_mutex_unlock(&cache.bounds_mtx);

            initialize_est_bounds(a2, b1, b2, bpar, &cb.beb);
            cb.beb.state = SET;
            
            pthread_mutex_lock(&cache.bounds_mtx);
            k = kh_get(bounds_hash, cache.bounds, bk.val);
            if (k == kh_end(cache.bounds))
                k = kh_put(bounds_hash, cache.bounds, bk.val, &ret);
                
            kh_val(cache.bounds, k) = cb;
            pthread_mutex_unlock(&cache.bounds_mtx);
            pthread_cond_signal(&cache.bounds_cond);
            *state = locate_cell(&cb.beb, a1);
            state_unset = 0;
        }
        else
        {
            /* some other thread has already inserted this key.  but,
               that thread may still be working on the value... */
            cb = kh_val(cache.bounds, k);
            while (cb.beb.state == PENDING)
            {
                pthread_cond_wait(&cache.bounds_cond, &cache.bounds_mtx);
                k = kh_get(bounds_hash, cache.bounds, bk.val);
                cb = kh_val(cache.bounds, k);
            }
            ++cb.hitcount;
            kh_val(cache.bounds, k) = cb;
            pthread_mutex_unlock(&cache.bounds_mtx);
            *state = locate_cell(&cb.beb, a1);
            state_unset = 0;
        }
    }
    if (! hash_frozen) winnow_bounds_hash(1);
    return state_unset;
}

    
/* test two dirichlets based on their counts. if the pattern of counts
   is cacheable, 'cacheable' is set to 1, and the cache is queried for
   the entry, and cache_was_set is set to 1 if found. */
enum fuzzy_state
cached_dirichlet_diff(unsigned *a,
                      unsigned *b,
                      struct binomial_est_params *bpar,
                      unsigned *cacheable,
                      unsigned *cache_was_set)
{
    unsigned lim[] = { cache.max_depth, cache.max_depth, 1, 1 };
    unsigned p[4];
    enum fuzzy_state state = UNCHANGED;
    int state_unset;

    find_cacheable_permutation(a, b, lim, p, cacheable);
    if (*cacheable)
    {
        unsigned a1 = a[p[0]], a2 = a[p[1]], b1 = b[p[0]], b2 = b[p[1]];
        state_unset = 
            get_fuzzy_state(bpar, a1, a2, b1, b2, 
                            &state, bpar->bounds_hash_frozen, cache_was_set);
    }
    else state_unset = 1;

    if (state_unset)
    {
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
