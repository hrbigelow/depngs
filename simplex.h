#ifndef _SIMPLEX_H
#define _SIMPLEX_H

/* The regular 4-simplex with unit-length sides, and centered about
   the origin, has corner points (let R = 1/(2 * sqrt(2))
   ( 1/2, 0, -R ), ( -1/2, 0, -R ), ( 0, 1/2, R ), ( 0, -1/2, R )
*/

double simplex_corners[4][3];


/* convert coordinates between the 4-simplex with unit side (which is
   embedded in R3) to barycentric coordinates (in R4, normalized,
   positive coordinates, weights on the 4 corners of the 4-simplex) */
void barycentric_to_simplex(double *b, double *s);
void simplex_to_barycentric(double *s, double *b);

/* compuates the square of the euclidean distance between a and b,
   points in R3 */
double square_dist3(double *a, double *b);

/* euclidean distance of two points in R3 */
double dist3(double *a, double *b);

/* return 1 if x is inside the 4-simplex */
unsigned inside_simplex(double *x);

#endif /* _SIMPLEX_H */
