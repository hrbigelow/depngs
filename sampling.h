#ifndef _SAMPLING_H
#define _SAMPLING_H

#include <stdlib.h>
#include <stdio.h>

/* Group of functions to generate individual samples from statistical
   distributions */
void compute_marginal_quantiles(double *sample_points,
                                size_t n_points,
                                size_t n_dims,
                                size_t sort_dimension,
                                const double *quantiles,
                                size_t n_quantiles,
                                double *quantile_values);

size_t marginal_quantiles_locus_bytes(size_t n_quantiles);

char *print_marginal_quantiles(char *out_buf,
                               double *sample_points,
                               size_t n_points,
                               const char *line_label,
                               const char *sums_label,
                               const double *quantiles, 
                               size_t n_quantiles);

void compute_marginal_wquantiles(double *sample_points,
                                 double *weights,
                                 size_t n_points,
                                 size_t n_dims,
                                 size_t sort_dimension,
                                 const double *quantiles,
                                 size_t n_quantiles,
                                 double *quantile_values);

char *print_marginal_wquantiles(char *out_buf, 
                                double *sample_points,
                                double *weights,
                                size_t n_points,
                                const char *line_label,
                                const double **quantile_values,
                                size_t n_quantiles);


void print_numerical_cdfs(FILE *out_fh, 
                          const char *label,
                          double *sample_points,
                          size_t n_points,
                          size_t dim);

#endif // _SAMPLING_H
