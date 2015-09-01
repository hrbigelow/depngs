#ifndef _DIR_CACHE_H
#define _DIR_CACHE_H

/* Functions for compiling a cache of both dirichlet point sets and
   binomial estimation bounds that are abundant in the data. */

#include "defs.h"
#include "khash.h"
#include "binomial_est.h"

struct dir_cache_params {
    unsigned n_bounds;
    unsigned min_ct_keep_bound;
    unsigned n_point_sets;
    unsigned max_sample_points;
    unsigned long n_max_survey_loci;
    const char *fasta_file; /* needed to initialize batch_pileup */
};

/* call once at start of program */
void
dir_cache_init(struct dir_cache_params dc_par);


/* call once at end of program */
void
dir_cache_free();


/* attempt to return the stored points from alpha counts, NULL if not
   found. */
POINT *
dir_cache_try_get_points(unsigned *alpha);

/* return a cached bounds or NULL if not cached */
struct binomial_est_bounds *
dir_cache_try_get_bounds(unsigned a2, unsigned b1, unsigned b2);


/* main routine for running the survey.  read chunks of input until
   accumulating at least n_bounds and n_point_sets.  The number of
   occurrences of each bound is recorded as well, and only those with
   >= min_ct_keep_bound are counted towards the n_bounds
   requirement. */
void
run_survey(struct bam_filter_params bf_par,
           struct bam_scanner_info *reader_buf,
           unsigned pseudo_depth,
           unsigned long n_max_survey_loci,
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
