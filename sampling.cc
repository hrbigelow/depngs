#include "sampling.h"
#include "defs.h"
#include <string.h>
#include <algorithm>
#include <math.h>
// #include <map>
// #include <cmath>

#ifdef __cplusplus
extern "C" {
#endif 

/* compute marginal quantiles on a given dimension of the sample
   points */
void
compute_marginal_quantiles(double *sample_points,
                           size_t n_points,
                           size_t n_dims,
                           size_t sort_dimension,
                           const double *quantiles,
                           size_t n_quantiles,
                           double *quantile_values)
{
    /* copy appropriate dimension */
    double *dim_points = (double *)malloc(sizeof(double) * n_points);
    double *start = dim_points, *end = start + n_points;
    double *cut;
    
    double *ps = sample_points + sort_dimension;
    double *p = dim_points;
    for ( ; p != end; ++p, ps += n_dims) *p = *ps;
    
    size_t f;
    for (f = 0; f != n_quantiles; ++f) {
        cut = dim_points + (size_t)round(quantiles[f] * n_points);
        std::nth_element(start, cut, end);
        quantile_values[f] = cut == end ? 0.0 : *cut;
        start = cut;
    }
    free(dim_points);
}


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


/* given a set of n_points weighted points and sum of those weights,
   partition the range into two regions with the coordinate in the
   first region always less than that in the second (using
   nth_element) and whose mass of region A < quantile * wgt_sum. */
struct weighted_coord *
weighted_quantiles_aux(struct weighted_coord *wgt_points,
                       size_t n_points,
                       double sum_wgt,
                       double quantile)
{
    struct weighted_coord 
        *wb = wgt_points, *we = wb + n_points, *cut;

    /* Find an initial estimate */
    struct less_wcoord lesswc;
    double target_wgt = quantile * sum_wgt;
    cut = wb + (size_t)floor(quantile * (double)n_points);
    std::nth_element(wb, cut, we, lesswc);

    if (cut == wb || cut == we) return cut;

    /* If cut is somewhere in the middle, then we have two non-empty
       intervals to possibly traverse.  */
    double sum_wgt_lft, sum_wgt_rgt;

    struct weighted_coord *wp;
    if (cut - wb < we - cut) {
        for (wp = wb, sum_wgt_lft = 0; wp != cut; ++wp)
            sum_wgt_lft += wp->weight;
        sum_wgt_rgt = sum_wgt - sum_wgt_lft;
    } else {
        for (wp = cut, sum_wgt_rgt = 0; wp != we; ++wp)
            sum_wgt_rgt += wp->weight;
        sum_wgt_lft = sum_wgt - sum_wgt_rgt;
    }
    /* choose to recurse on left or right sub-interval */
    if (sum_wgt_lft < target_wgt) {
        n_points = we - cut;
        wb = cut;
        quantile = (target_wgt - sum_wgt_lft) / sum_wgt_rgt;
        sum_wgt = sum_wgt_rgt;
    }
    else {
        n_points = cut - wb;
        quantile = target_wgt / sum_wgt_lft;
        sum_wgt = sum_wgt_lft;
    }
    return weighted_quantiles_aux(wb, n_points, sum_wgt, quantile);
}


/* computes the marginal quantiles from a set of weighted sample
   points, selecting the 'dim' component of the point. */
void
compute_marginal_wquantiles(double *sample_points,
                            double *weights,
                            size_t n_points,
                            size_t n_dims,
                            size_t dim,
                            const double *quantiles,
                            size_t n_quantiles,
                            double *quantile_values)
{
    /* copy appropriate dimension */
    struct weighted_coord 
        *wgt_points = 
        (struct weighted_coord *)malloc(sizeof(struct weighted_coord) * n_points);
    
    struct weighted_coord *start = wgt_points, *end = start + n_points, *cut;

    double *ps = sample_points + dim;
    double *pw = weights;
    double sum_wgt = 0;
    for (cut = start; cut != end; ++cut, ps += n_dims, ++pw) {
        cut->coord = *ps;
        cut->weight = *pw;
        sum_wgt += *pw;
    }

    size_t f;
    for (f = 0; f != n_quantiles; ++f) {
        cut = weighted_quantiles_aux(wgt_points, n_points, sum_wgt, quantiles[f]);
        quantile_values[f] = cut->coord;
    }

    free(wgt_points);
}


double compute_marginal_mean(double *points,
                             double *weights,
                             size_t n_points,
                             size_t n_dims,
                             size_t dim)
{
    double 
        *pb,
        *wb = weights, *we = wb + n_points, 
        sum_wgts = 0, marg = 0;

    while (wb != we) sum_wgts += *wb++;

    for (wb = weights, pb = points + dim; wb != we; ++wb, pb += n_dims) 
        marg += *pb * *wb;

    marg /= sum_wgts;
    return marg;
}
                             

#if 0
/* prints out return next write position after writing to out_buf */
char *print_marginal_quantiles(char *out_buf, 
                               const char *line_label,
                               const double **quantile_values,
                               const double *means,
                               size_t n_quantiles)
{
    std::multimap<double, size_t, std::greater<double> > dim_to_mean;
    unsigned d, q;
    for (d = 0; d != NUM_NUCS; ++d)
        dim_to_mean.insert(std::make_pair(means[d], d));

    /* calculate mean rank order */
    size_t mean_rank_order[NUM_NUCS];
    std::multimap<double, size_t, std::greater<double> >::iterator dtm;
    d = 0;
    for (dtm = dim_to_mean.begin(); dtm != dim_to_mean.end(); ++dtm)
        mean_rank_order[(*dtm).second] = d++;

    static const char dimension_labels[] = "ACGT";
    for (d = 0; d != NUM_NUCS; ++d)
    {
        out_buf += 
            sprintf(out_buf, "%s\t%c\t%Zu\t%10.8f", line_label, 
                    dimension_labels[d], mean_rank_order[d], means[d]);

        for (q = 0; q != n_quantiles; ++q)
            out_buf += sprintf(out_buf, "\t%10.8f", quantile_values[d][q]);

        out_buf += sprintf(out_buf, "\n");
    }
    
    return out_buf;
}
#endif


#if 0
// return next write position after writing to out_buf
char *print_marginal_quantiles(char *out_buf, 
                               double *sample_points,
                               size_t n_points,
                               const char *line_label,
                               const char *sums_label,
                               const double *quantiles, 
                               size_t n_quantiles)
{
    double quantile_sums[MAX_NUM_QUANTILES], quantile_values[MAX_NUM_QUANTILES];
    memset(quantile_sums, 0, sizeof(quantile_sums));

    double mean[] = { 0, 0, 0, 0 }, mean_sum = 0.0;

    /* calculate mean */
    double *point = sample_points;
    size_t d, q, p = 0;
    for ( ; p != n_points; ++p, point += NUM_NUCS)
        for (d = 0; d != NUM_NUCS; ++d) mean[d] += point[d];

    double n_points_inv = 1.0 / (double)n_points;
    std::multimap<double, size_t, std::greater<double> > dim_to_mean;
    for (d = 0; d != NUM_NUCS; ++d)
    {
        mean[d] *= n_points_inv;
        dim_to_mean.insert(std::make_pair(mean[d], d));
        mean_sum += mean[d];
    }

    /* calculate mean rank order */
    size_t mean_rank_order[NUM_NUCS];
    std::multimap<double, size_t, std::greater<double> >::iterator dtm;
    d = 0;
    for (dtm = dim_to_mean.begin(); dtm != dim_to_mean.end(); ++dtm)
        mean_rank_order[(*dtm).second] = d++;

    static const char dimension_labels[] = "ACGT";
    for (d = 0; d != NUM_NUCS; ++d)
    {
        out_buf += 
            sprintf(out_buf, "%s\t%c\t%Zu\t%10.8f", line_label, 
                    dimension_labels[d], mean_rank_order[d], mean[d]);

        compute_marginal_quantiles(sample_points, n_points, NUM_NUCS, d, 
                                   quantiles, n_quantiles, quantile_values);
        
        for (q = 0; q != n_quantiles; ++q)
        {
            out_buf += sprintf(out_buf, "\t%10.8f", quantile_values[q]);
            quantile_sums[q] += quantile_values[q];
        }
        out_buf += sprintf(out_buf, "\n");
    }
    
    out_buf += sprintf(out_buf, "%s\t%s\t%s\t%10.8f", 
                       line_label, sums_label, sums_label, mean_sum);
    
    for (q = 0; q != n_quantiles; ++q)
        out_buf += sprintf(out_buf, "\t%10.8f", quantile_sums[q]);
    
    out_buf += sprintf(out_buf, "\n");
    return out_buf;
}
#endif



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
                          size_t n_points,
                          size_t ndim)
{
    // reshape sample_points so each dimension is held together
    double **dim_points = new double*[n_points * ndim];
    size_t *ranks = new size_t[n_points * ndim];
    size_t d;

    double **p = dim_points;
    size_t *r = ranks;

    for (d = 0; d != ndim; ++d)
    {
        double *ps = sample_points + d;
        double **end = p + n_points;
        for ( ; p != end; ++p, ps += 4)
            *p = ps;
    }
    
    p = dim_points;
    less_ptr lp;
    for (d = 0; d != ndim; ++d, p += n_points)
        std::sort(p, p + n_points, lp);

    // compute ranks
    p = dim_points;
    r = ranks;
    for (d = 0; d != ndim; ++d, r += n_points)
        for (size_t pi = 0; pi != n_points; ++pi)
            r[(*p++ - (sample_points + d)) / 4] = pi;
    
    for (size_t pi = 0; pi != n_points; ++pi)
    {
        fprintf(out_fh, "%Zu\t%s", pi, label);
        for (d = 0; d != ndim; ++d)
            fprintf(out_fh, "\t%20.18f", sample_points[ndim * pi + d]);

        for (d = 0; d != ndim; ++d)
            fprintf(out_fh, "\t%Zu", ranks[n_points * d + pi]);

        fprintf(out_fh, "\n");
    }
    delete dim_points;
    delete ranks;
}


#ifdef __cplusplus
}
#endif 
