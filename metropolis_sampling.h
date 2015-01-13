#ifndef _METROPOLIS_SAMPLING_H
#define _METROPOLIS_SAMPLING_H

#include "nucleotide_stats.h"

void tune_proposal(struct packed_counts *cts,
                   struct posterior_settings *set,
                   const double *prior_alpha,
                   double *proposal_alpha,
                   double *estimated_mean,
                   double *points_buf);

void metropolis_sampling(unsigned short start_point, unsigned short n_points_wanted,
                         const struct packed_counts *cts,
                         const double *logu, /* logs of U[0, 1] values */
                         double *proposal_alpha,
                         size_t nth, /* collect every nth point */
                         double *sample_points);

#endif /* METROPOLIS_SAMPLING_H */
