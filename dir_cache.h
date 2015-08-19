#ifndef _DIR_CACHE_H
#define _DIR_CACHE_H

#include "defs.h"
#include "khash.h"
#include "binomial_est.h"

/* use union alpha_large_key as key */
KHASH_MAP_INIT_INT64(points_h, POINT *);

/* use union bounds_key as key */
KHASH_MAP_INIT_INT64(bounds_h, struct binomial_est_bounds);

khash_t(points_h) *g_points_hash;
khash_t(bounds_h) *g_bounds_hash;

                  // {};

struct dir_cache_params {
    unsigned n_bounds;
    unsigned min_ct_keep_bound;
    unsigned n_point_sets;
    unsigned max_sample_points;
    const char *fasta_file; /* needed to initialize batch_pileup */
};

/* call once at start of program */
void
dir_cache_init(struct dir_cache_params dc_par);


/* call once at end of program */
void
dir_cache_free();


/* main routine for running the survey.  read chunks of input until
   accumulating at least n_bounds and n_point_sets.  The number of
   occurrences of each bound is recorded as well, and only those with
   >= min_ct_keep_bound are counted towards the n_bounds
   requirement. */
void
run_survey(void **reader_pars,
           struct contig_region *qbeg,
           struct contig_region *qend,
           unsigned n_threads,
           unsigned n_max_reading,
           unsigned long max_input_mem);


/* populates the internal hash g_points_hash (and backing buffer
   g_point_sets_buf) with dirichlet points parameterized by the keys
   in g_pt_hash.  run_survey and dirichlet_points_gen_init to be
   called first. */
void
generate_point_sets(unsigned n_threads);


/* populate g_bounds_hash with computed est bounds for the bounds
   tuples surveyed. */
void
generate_est_bounds(unsigned n_threads);


#endif /* _DIR_CACHE_H */
