#include "sampling.h"
#include "defs.h"
#include <string.h>
#include <algorithm>
#include <map>
#include <cmath>

/* compute marginal quantiles on a given dimension of the sample
   points */
void compute_marginal_quantiles(double *sample_points,
                                size_t num_points,
                                size_t sort_dimension,
                                const double *quantiles,
                                size_t num_quantiles,
                                double *quantile_values)
{
    /* copy appropriate dimension */
    double *dim_points = (double *)malloc(sizeof(double) * num_points);
    double *start = dim_points, *end = start + num_points;
    double *cut;

    double *ps = sample_points + sort_dimension;
    double *p = dim_points;
    for ( ; p != end; ++p, ps += 4) *p = *ps;

    size_t f;
    for (f = 0; f != num_quantiles; ++f)
    {
        cut = dim_points + (size_t)std::round(quantiles[f] * num_points);
        std::nth_element(start, cut, end);
        quantile_values[f] = cut == end ? 0.0 : *cut;
        start = cut;
    }
    free(dim_points);
}

#define NDIM 4


struct weighted_coord {
    double coord;
    double weight;
};

struct less_wcoord {
    bool operator()(const struct weighted_coord &c1,
                    const struct weighted_coord &c2)
    {
        return c1.coord < c2.coord;
    }
};

void compute_marginal_wquantiles(double *sample_points,
                                 double *weights,
                                 size_t n_points,
                                 size_t n_dims,
                                 size_t sort_dimension,
                                 const double *quantiles,
                                 size_t n_quantiles,
                                 double *quantile_values)
{
    /* copy appropriate dimension */
    struct weighted_coord 
        *wgt_points = 
        (struct weighted_coord *)malloc(sizeof(struct weighted_coord) * n_points);
    
    struct weighted_coord *start = wgt_points, *end = start + n_points, *cut;

    double *ps = sample_points + sort_dimension;
    double *pw = weights;
    for (cut = start; cut != end; ++cut, ps += n_dims, ++pw) 
    {
        cut->coord = *ps;
        cut->weight = *pw;
    }

    size_t f;
    double running_sum = 0;
    struct less_wcoord lesswc;
    for (f = 0; f != n_quantiles; ++f)
    {
        cut = wgt_points + (size_t)std::round(quantiles[f] * n_points);
        std::nth_element(start, cut, end, lesswc);
        while (start != cut) running_sum += start++->weight;
        quantile_values[f] = running_sum;
    }

    /* normalize */
    double running_sum_inv = 1.0 / running_sum;
    for (f = 0; f != n_quantiles; ++f)
        quantile_values[f] *= running_sum_inv;

    free(wgt_points);
}


// return next write position after writing to out_buf
char *print_marginal_wquantiles(char *out_buf, 
                                double *sample_points,
                                double *weights,
                                size_t num_points,
                                const char *line_label,
                                const double **quantile_values,
                                size_t num_quantiles)
{
    double mean[] = { 0, 0, 0, 0 };

    /* calculate mean */
    double *point = sample_points;
    double weights_sum = 0;
    size_t d, q, p = 0;
    for ( ; p != num_points; ++p, point += NDIM)
    {
        weights_sum += weights[p];
        for (d = 0; d != NDIM; ++d) mean[d] += point[d] * weights[p];
    }
    double weights_sum_inv = 1.0 / weights_sum;

    std::multimap<double, size_t, std::greater<double> > dim_to_mean;
    for (d = 0; d != NDIM; ++d)
    {
        mean[d] *= weights_sum_inv;
        dim_to_mean.insert(std::make_pair(mean[d], d));
    }

    /* calculate mean rank order */
    size_t mean_rank_order[NDIM];
    std::multimap<double, size_t, std::greater<double> >::iterator dtm;
    d = 0;
    for (dtm = dim_to_mean.begin(); dtm != dim_to_mean.end(); ++dtm)
        mean_rank_order[(*dtm).second] = d++;

    static const char dimension_labels[] = "ACGT";
    for (d = 0; d != NDIM; ++d)
    {
        out_buf += 
            sprintf(out_buf, "%s\t%c\t%Zu\t%10.8f", line_label, 
                    dimension_labels[d], mean_rank_order[d], mean[d]);

        for (q = 0; q != num_quantiles; ++q)
            out_buf += sprintf(out_buf, "\t%10.8f", quantile_values[d][q]);

        out_buf += sprintf(out_buf, "\n");
    }
    
    return out_buf;
}


// return next write position after writing to out_buf
char *print_marginal_quantiles(char *out_buf, 
                               double *sample_points,
                               size_t num_points,
                               const char *line_label,
                               const char *sums_label,
                               const double *quantiles, 
                               size_t num_quantiles)
{
    double quantile_sums[MAX_NUM_QUANTILES], quantile_values[MAX_NUM_QUANTILES];
    memset(quantile_sums, 0, sizeof(quantile_sums));

    double mean[] = { 0, 0, 0, 0 }, mean_sum = 0.0;

    /* calculate mean */
    double *point = sample_points;
    size_t d, q, p = 0;
    for ( ; p != num_points; ++p, point += NDIM)
        for (d = 0; d != NDIM; ++d) mean[d] += point[d];

    double num_points_inv = 1.0 / (double)num_points;
    std::multimap<double, size_t, std::greater<double> > dim_to_mean;
    for (d = 0; d != NDIM; ++d)
    {
        mean[d] *= num_points_inv;
        dim_to_mean.insert(std::make_pair(mean[d], d));
        mean_sum += mean[d];
    }

    /* calculate mean rank order */
    size_t mean_rank_order[NDIM];
    std::multimap<double, size_t, std::greater<double> >::iterator dtm;
    d = 0;
    for (dtm = dim_to_mean.begin(); dtm != dim_to_mean.end(); ++dtm)
        mean_rank_order[(*dtm).second] = d++;

    static const char dimension_labels[] = "ACGT";
    for (d = 0; d != NDIM; ++d)
    {
        out_buf += 
            sprintf(out_buf, "%s\t%c\t%Zu\t%10.8f", line_label, 
                    dimension_labels[d], mean_rank_order[d], mean[d]);

        compute_marginal_quantiles(sample_points, num_points, d, 
                                   quantiles, num_quantiles, quantile_values);
        
        for (q = 0; q != num_quantiles; ++q)
        {
            out_buf += sprintf(out_buf, "\t%10.8f", quantile_values[q]);
            quantile_sums[q] += quantile_values[q];
        }
        out_buf += sprintf(out_buf, "\n");
    }
    
    out_buf += sprintf(out_buf, "%s\t%s\t%s\t%10.8f", 
                       line_label, sums_label, sums_label, mean_sum);
    
    for (q = 0; q != num_quantiles; ++q)
        out_buf += sprintf(out_buf, "\t%10.8f", quantile_sums[q]);
    
    out_buf += sprintf(out_buf, "\n");
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
