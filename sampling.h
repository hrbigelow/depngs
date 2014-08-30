#ifndef _SAMPLING_H
#define _SAMPLING_H

#include "tools.h"

/*
  Group of functions to generate individual samples from statistical distributions
 */


void compute_marginal_quantiles(double * sample_points,
                                size_t num_points,
                                size_t sort_dimension,
                                double const* quantiles,
                                size_t num_quantiles,
                                double * quantile_values);


size_t marginal_quantiles_locus_bytes(size_t num_quantiles);

char * print_marginal_quantiles(char * out_buf,
                                double * sample_points,
                                size_t num_points,
                                double const* mode_point,
                                char const* line_label,
                                char const** dimension_labels,
                                char const* sums_label,
                                double const* quantiles, 
                                size_t num_quantiles);

void print_numerical_cdfs(FILE * out_fh, 
                          char const* label,
                          double * sample_points,
                          size_t num_points,
                          size_t dim);

#endif // _SAMPLING_H
