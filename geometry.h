#ifndef _GEOMETRY_H
#define _GEOMETRY_H



/* computes the euclidean distance between corresponding pairs of
   points with n_dims, populating square_dist_buf */
void compute_square_dist(const double *points1,
                         const double *points2,
                         unsigned n_points,
                         unsigned n_dims,
                         double *square_dist_buf);

#endif /* _GEOMETRY_H */
