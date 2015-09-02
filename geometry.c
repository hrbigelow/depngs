#include "geometry.h"
#include "yepCore.h"

#include <stdlib.h>

/* computes the euclidean distance between corresponding pairs of
   points with n_dims, populating square_dist_buf */
void compute_square_dist(const double *points1,
                         const double *points2,
                         unsigned n_points,
                         unsigned n_dims,
                         double *square_dist_buf)
{
    unsigned n_comp = n_points * n_dims;
    double *diff = (double *)malloc(n_comp * sizeof(double));
    double *diffsq = (double *)malloc(n_comp * sizeof(double));
    (void)yepCore_Subtract_V64fV64f_V64f(points1, points2, diff, n_comp);
    (void)yepCore_Multiply_V64fV64f_V64f(diff, diff, diffsq, n_comp);

    unsigned p, c;
    double *pt = diffsq;
    for (p = 0; p != n_points; ++p) {
        square_dist_buf[p] = 0;
        for (c = 0; c != n_dims; ++c)
            square_dist_buf[p] += *pt++;
    }
    free(diff);
    free(diffsq);
}
