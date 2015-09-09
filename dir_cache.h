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
    unsigned pseudo_depth;
    double prior_alpha;
    const char *fasta_file; /* needed to initialize batch_pileup */
};

/* call once at start of program */
void
dir_cache_init(struct binomial_est_params be_par,
               struct dir_cache_params dc_par,
               struct bam_filter_params bf_par,
               struct bam_scanner_info *reader_buf,
               unsigned n_max_reading,
               unsigned long max_input_mem,
               unsigned n_threads);


/* call once at end of program */
void
dir_cache_free();


/* attempt to return the stored points from alpha counts, NULL if not
   found. */
POINT *
dir_cache_try_get_points(unsigned *alpha);


struct dir_cache_pair_key {
    uint64_t key;
    unsigned char is_valid;
    unsigned char is_ref_change;
};


/* try to get the fuzzy_state corresponding to this key, or return
   STATE_UNKNOWN. */
enum fuzzy_state
dir_cache_try_get_state(struct dir_cache_pair_key key);



/* classify a pair of alphas with their permutation provided */
struct dir_cache_pair_key
dir_cache_classify_alphas(const unsigned *alpha1,
                          const unsigned *alpha2);


/* provide a fuzzy_state change by calculation.  called when the cache
   doesn't have an entry. */
enum fuzzy_state
dir_cache_calc_state(const unsigned *alpha1,
                     const unsigned *alpha2,
                     struct dir_points *dist1,
                     struct dir_points *dist2);


/* return index of pseudo_alpha component, or -1 if this is not a
   pseudo_alpha. */
int
dir_cache_is_pseudo_alpha(const unsigned *alpha);




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
/* void */
/* generate_est_bounds(unsigned n_threads); */


/* populate g_ref_change_hash with computed bounds for sample-to-REF
   comparisons. */
void
generate_ref_change(unsigned n_threads);


/* populate g_sam_change_hash with computed changes for
   sample-to-sample comparisons. */
void
generate_sam_change(unsigned n_threads);

#endif /* _DIR_CACHE_H */
