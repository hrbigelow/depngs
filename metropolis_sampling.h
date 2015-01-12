#ifndef _METROPOLIS_SAMPLING_H
#define _METROPOLIS_SAMPLING_H

void metropolis_sampling(unsigned short start_point, unsigned short n_points_wanted,
                         const struct cpd_count *term, size_t n_term,
                         double *proposal_alpha,
                         double *prior_alpha,
                         double *logu, /* logs of U[0, 1] values */
                         size_t nth, /* collect every nth point */
                         double *sample_points);

#endif /* METROPOLIS_SAMPLING_H */
