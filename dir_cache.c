/* Utilities for caching bounds and points */

/* use union alpha_large_key as key */
KHASH_MAP_INIT_INT64(points_tuple_h, unsigned);

/* use union bounds_key as key */
KHASH_MAP_INIT_INT64(bounds_tuple_h, unsigned);

/* use union alpha_large_key as key */
KHASH_MAP_INIT_INT64(points_h, double *);

/* use union bounds_key as key */
KHASH_MAP_INIT_INT64(bounds_h, struct binomial_est_bounds);


/* Survey */
static khash_t(points_tuple_h) *g_points_tuple_hash;
static khash_t(bounds_tuple_h) *g_bounds_tuple_hash;

static __thread khash_t(points_tuple_h) *t_points_tuple_hash;
static __thread khash_t(bounds_tuple_h) *t_bounds_tuple_hash;

/* Prepopulation */
static khash_t(points_h) *g_points_hash;
static khash_t(bounds_h) *g_bounds_hash;

static __thread khash_t(points_h) *t_points_hash;
static __thread khash_t(bounds_h) *t_bounds_hash;

/* Operation */



