/* compute quantiles from weighted points */

#include "wquantile.h"
#include "ksort.h"

#include <math.h>

struct weighted_coord {
    double coord;
    double weight;
};

static int
less_coord(struct weighted_coord a, struct weighted_coord b)
{
    return a.coord < b.coord;
}


KSORT_INIT(wpoint_sort, struct weighted_coord, less_coord);


void
compute_marginal_wgt_quantiles(double *sample_points,
                               double *weights,
                               unsigned n_points,
                               unsigned n_dims,
                               const double *quantiles,
                               unsigned n_quantiles,
                               quantile_vals_t *quantile_values)
{
    struct weighted_coord
        *wcb = malloc(sizeof(struct weighted_coord) * n_points),
        *wce = wcb + n_points,
        *wc;

    double *ps = sample_points;
    double *pw = weights, *pe = pw + n_points;
    double sum_wgt = 0;
    for (pw = weights; pw != pe; ++pw)
        sum_wgt += *pw;

    quantile_vals_t quantiles_adj;
    unsigned q;
    for (q = 0; q != n_quantiles; ++q)
        quantiles_adj[q] = quantiles[q] * sum_wgt;
    
    unsigned d;
    for (d = 0; d != n_dims; ++d) {
        /* collate the sort dimension with the weights */
        for (wc = wcb, ps = sample_points + d, pw = weights;
             wc != wce;
             ++wc, ps += n_dims, ++pw)
            *wc = (struct weighted_coord){ *ps, *pw };

        ks_introsort(wpoint_sort, n_points, wcb);
        double cumul_wgt = 0;
        
        wc = wcb;
        for (q = 0; q != n_quantiles; ++q) {
            while (cumul_wgt < quantiles_adj[q] && wc != wce)
                cumul_wgt += wc++->weight;
           /* use coordinate of last point if quantile of 1 is requested */
            quantile_values[d][q] = wc == wce ? (wc - 1)->coord : wc->coord;
        }
    }
    free(wcb);
}


/* Compute the requested set of distance quantile values from two sets
   of weighted points. */
void
compute_dist_wgt_quantiles(double *square_dist_buf,
                           double *weights_buf,
                           unsigned n_points,
                           const double *quantiles,
                           unsigned n_quantiles,
                           quantile_vals_t dist_quantile_values)
{
    compute_marginal_wgt_quantiles(square_dist_buf, weights_buf, n_points, 1,
                                   quantiles, n_quantiles,
                                   (quantile_vals_t *)dist_quantile_values);

    unsigned q;
    for (q = 0; q != n_quantiles; ++q) 
        dist_quantile_values[q] = sqrt(dist_quantile_values[q]);
}
