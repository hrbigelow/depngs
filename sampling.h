#ifndef _SAMPLING_H
#define _SAMPLING_H

#include <stdlib.h>
#include <stdio.h>

/* Group of functions to generate individual samples from statistical
   distributions */

#ifdef __cplusplus
extern "C" {
#endif

void
compute_marginal_quantiles(double *sample_points,
                           size_t n_points,
                           size_t n_dims,
                           size_t sort_dimension,
                           const double *quantiles,
                           size_t n_quantiles,
                           double *quantile_values);


void
compute_marginal_wquantiles(double *sample_points,
                            double *weights,
                            size_t n_points,
                            size_t n_dims,
                            size_t sort_dimension,
                            const double *quantiles,
                            size_t n_quantiles,
                            double *quantile_values);

double compute_marginal_mean(double *points,
                             double *weights,
                             size_t n_points,
                             size_t n_dims,
                             size_t dim);


char *print_marginal_quantiles(char *out_buf, 
                               const char *line_label,
                               const double **quantile_values,
                               const double *means,
                               size_t n_quantiles);


void print_numerical_cdfs(FILE *out_fh, 
                          const char *label,
                          double *sample_points,
                          size_t n_points,
                          size_t dim);

#ifdef __cplusplus
}
#endif

#endif // _SAMPLING_H
