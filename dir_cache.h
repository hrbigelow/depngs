#ifndef _DIR_CACHE_H
#define _DIR_CACHE_H

/* populates the internal hash g_points_hash (and backing buffer
   g_point_sets_buf) with dirichlet points parameterized by the keys
   in g_pt_hash.  run_survey and dirichlet_points_gen_init to be
   called first. */
void
generate_point_sets(unsigned n_threads,
                    unsigned max_sample_points);




#endif /* _DIR_CACHE_H */
