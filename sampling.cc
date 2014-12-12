#include "sampling.h"

#include <algorithm>

// compute marginal quantiles on a given dimension of the sample
// points
void compute_marginal_quantiles(double *sample_points,
                                size_t num_points,
                                size_t sort_dimension,
                                const double *quantiles,
                                size_t num_quantiles,
                                double *quantile_values)
{
    // copy appropriate dimension
    double *dim_points = new double[num_points];
    double *start = dim_points, *end = start + num_points;
    double *cut;

    double *ps = sample_points + sort_dimension;
    double *p = dim_points;
    for ( ; p != end; ++p, ps += 4)
    {
        *p = *ps;
    }

    for (size_t f = 0; f != num_quantiles; ++f)
    {
        cut = dim_points + static_cast<size_t>(std::round(quantiles[f] *num_points));
        std::nth_element(start, cut, end);
        quantile_values[f] = cut == end ? 0.0 : *cut;
        start = cut;
    }
    delete dim_points;

}

// this should be in synch with print_marginal_quantiles
// there are 5 lines per locus
size_t marginal_quantiles_locus_bytes(size_t num_quantiles)
{
    return 5 * (71 + (10 * num_quantiles) + 11 + num_quantiles);
}

// return next write position after writing to out_buf
char *print_marginal_quantiles(char *out_buf, 
                               double *sample_points,
                               size_t num_points,
                               const char *line_label,
                               const char **dimension_labels,
                               const char *sums_label,
                               const double *quantiles, 
                               size_t num_quantiles)
{

    size_t num_dimensions = 4;

    double *quantile_sums = new double[num_quantiles];
    std::fill(quantile_sums, quantile_sums + num_quantiles, 0.0);

    double *quantile_values = new double[num_quantiles];

    double mean_sum = 0.0;

    double *mean = new double[num_dimensions];
    std::multimap<double, size_t, std::greater<double> > dim_to_mean;

    //calculate mean
    double *point = sample_points;
    size_t p = 0;
    std::fill(mean, mean + num_dimensions, 0.0);
    for ( ; p != num_points; ++p, point += num_dimensions)
        for (size_t d = 0; d != num_dimensions; ++d)
            mean[d] += point[d];

    for (size_t d = 0; d != num_dimensions; ++d)
    {
        mean[d] /= num_points;
        dim_to_mean.insert(std::make_pair(mean[d], d));
        mean_sum += mean[d];
    }

    //calculate mean rank order
    size_t *mean_rank_order = new size_t[num_dimensions];
    std::multimap<double, size_t, std::greater<double> >::iterator dtm_iter;
    
    size_t d = 0;
    for (dtm_iter = dim_to_mean.begin(); dtm_iter != dim_to_mean.end(); 
         ++dtm_iter)
    {
        mean_rank_order[(*dtm_iter).second] = d;
        ++d;
    }

    for (size_t d = 0; d != num_dimensions; ++d)
    {
        out_buf += sprintf(out_buf, "%s\t%s\t%Zu", line_label, 
                           dimension_labels[d], mean_rank_order[d]);

        compute_marginal_quantiles(sample_points, num_points, d, 
                                   quantiles, num_quantiles, quantile_values);

        out_buf += sprintf(out_buf, "\t%10.8f", mean[d]);

        for (size_t q = 0; q != num_quantiles; ++q)
        {
            out_buf += sprintf(out_buf, "\t%10.8f", quantile_values[q]);
            quantile_sums[q] += quantile_values[q];
        }
        out_buf += sprintf(out_buf, "\n");
    }

    out_buf += sprintf(out_buf, "%s\t%s\t%s", line_label, sums_label, sums_label);
    out_buf += sprintf(out_buf, "\t%10.8f", mean_sum);

    for (size_t q = 0; q != num_quantiles; ++q)
    {
        out_buf += sprintf(out_buf, "\t%10.8f", quantile_sums[q]);
    }
    out_buf += sprintf(out_buf, "\n");
    
    delete mean_rank_order;
    delete mean;
    delete quantile_sums;
    delete quantile_values;

    return out_buf;
}

struct less_ptr
{
    bool operator()(double *a, double *b)
    {
        return (*a) < (*b);
    }
};


// transform a list of unordered tuples of (a,c,g,t) composition,
// and output their component-wise ranks as:
// i a_i, c_i, g_i, t_i, ra_i, rc_i, rg_i, rt_i
void print_numerical_cdfs(FILE *out_fh, 
                          const char *label,
                          double *sample_points,
                          size_t num_points,
                          size_t ndim)
{
    // reshape sample_points so each dimension is held together
    double **dim_points = new double*[num_points * ndim];
    size_t *ranks = new size_t[num_points * ndim];
    size_t d;

    double **p = dim_points;
    size_t *r = ranks;

    for (d = 0; d != ndim; ++d)
    {
        double *ps = sample_points + d;
        double **end = p + num_points;
        for ( ; p != end; ++p, ps += 4)
            *p = ps;
    }
    
    p = dim_points;
    less_ptr lp;
    for (d = 0; d != ndim; ++d, p += num_points)
        std::sort(p, p + num_points, lp);

    // compute ranks
    p = dim_points;
    r = ranks;
    for (d = 0; d != ndim; ++d, r += num_points)
        for (size_t pi = 0; pi != num_points; ++pi)
            r[(*p++ - (sample_points + d)) / 4] = pi;
    
    for (size_t pi = 0; pi != num_points; ++pi)
    {
        fprintf(out_fh, "%Zu\t%s", pi, label);
        for (d = 0; d != ndim; ++d)
            fprintf(out_fh, "\t%20.18f", sample_points[ndim * pi + d]);

        for (d = 0; d != ndim; ++d)
            fprintf(out_fh, "\t%Zu", ranks[num_points * d + pi]);

        fprintf(out_fh, "\n");
    }
    delete dim_points;
    delete ranks;
}
