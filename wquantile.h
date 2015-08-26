#ifndef _WQUANTILE_H
#define _WQUANTILE_H

#include "defs.h"

typedef double quantile_vals_t[MAX_NUM_QUANTILES];
typedef quantile_vals_t comp_quantile_vals_t[NUM_NUCS];

/* quantile_values[d][q] for d < n_dims and q < n_quantiles. */
void
compute_marginal_wgt_quantiles(double *sample_points,
                               double *weights,
                               unsigned n_points,
                               unsigned n_dims,
                               const double *quantiles,
                               unsigned n_quantiles,
                               quantile_vals_t *quantile_values);


/* wrapper for the 1-dimensional marginal weighted quantiles, that
   also transforms to square root. */
void
compute_dist_wgt_quantiles(double *square_dist_buf,
                           double *weights_buf,
                           unsigned n_points,
                           const double *quantiles,
                           unsigned n_quantiles,
                           quantile_vals_t dist_quantile_values);



#endif /* _WQUANTILE_H */
