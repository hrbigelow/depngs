#include "simplex.h"

#include <gsl/gsl_math.h>
#include <string.h>
#include <stdlib.h>
#include <stdio.h>

#include "yepCore.h"

/* 
   As it turns out, the determinant of the template is sqrt(2) /
   2. all we need to know is that it is positive.

   The Moore-Penrose Pseudoinverse equation holds that:
   
   Ax = b
   x = A'b + [I - A'A]w

   It turns out that the quantity [I - A'A] is the matrix with 0.25 in
   all cells.  Therefore, [I - A'A]w would equal the vector [0.25, 0.25,
   0.25, 0.25] regardless of the weight vector w (assuming it is
   normalized and positive)
*/

#define H (1.0 / 2.0)
#define R (1.0 / (2.0 * 1.41421356237309504880))
#define S (1.0 / 1.41421356237309504880)

static double simplex[][4] = { 
    { +H, -H, +0, +0 }, 
    { +0, +0, +H, -H }, 
    { -R, -R, +R, +R }
};

static double simpinv[][3] = { 
    { +1, +0, -S }, 
    { -1, +0, -S }, 
    { +0, +1, +S },  
    { +0, -1, +S } 
};

double simplex_corners[][3] = { 
    { +H, +0, -R }, 
    { -H, +0, -R }, 
    { +0, +H, +R }, 
    { +0, -H, +R } 
};

/* see http://steve.hollasch.net/cgindex/geometry/ptintet.html.  It
   outlines a technique for calculating whether a given point is
   within a tetrahedron.  The technique uses five determinants.  The
   first one is shown below.  The other four are constructed from this
   one by replacing the first three elements of one of the four rows
   with the coordinates of the point being tested.  Note: The
   determinant of this template is sqrt(2) / 2. But, all that is
   needed is to know that it is positive. See 'inside_simplex'
   below. */
static double det_template[][4] = {
    { +H, +0, -R, +1 },
    { -H, +0, -R, +1 },
    { +0, +H, +R, +1 },
    { +0, -H, +R, +1 }
};

void barycentric_to_simplex(double *b, double *s)
{
    unsigned r, c;
    for (r = 0; r != 3; ++r)
    {
        s[r] = 0;
        for (c = 0; c != 4; ++c)
            s[r] += simplex[r][c] * b[c];
    }
}


void simplex_to_barycentric(double *s, double *b)
{
    unsigned r, c;
    for (r = 0; r != 4; ++r)
    {
        b[r] = 0.25;
        for (c = 0; c != 3; ++c)
            b[r] += simpinv[r][c] * s[c];
    }
}

/* initialize the simplex matrix rows into repeated arrays */
static struct {
    double *row[3];
    unsigned n_rep;
} simplex_batch;

void init_simplex(unsigned n_rep)
{
    unsigned p, r;
    for (r = 0; r != 3; ++r)
    {
        simplex_batch.row[r] = malloc(n_rep * 4 * sizeof(double));
        for (p = 0; p != n_rep * 4; p += 4)
            memcpy(simplex_batch.row[r] + p, simplex[r], sizeof(simplex[r]));
    }
}


void free_simplex()
{
    unsigned r;
    for (r = 0; r != 3; ++r)
        free(simplex_batch.row[r]);
}

#if 0
/* calculate a batch of distances from a pair of barycentric
   coordinate arrays. break up the rows into 3 repeated arrays. */
void batch_bary_to_simplex_sqdist(const double *bary1, const double *bary2,
                                  double *sqdist, unsigned n_points)
{
    if(n_points > simplex_batch.n_rep)
    {
        fprintf(stderr, "Error: %s: number of points %u exceeds prepopulated"
                " number %u specified in 'init_simplex'.\n",
                __func__, n_points, simplex_batch.n_rep);
        exit(1);
    }
    unsigned n_comp = n_points * 3;
    double *prod = malloc(n_comp * sizeof(double));
    double *simp1 = calloc(n_comp, sizeof(double));
    double *simp2 = calloc(n_comp, sizeof(double));
    double *diff = malloc(n_comp * sizeof(double));
    double *diffsq = malloc(n_comp * sizeof(double));
    
    unsigned r;
    for (r = 0; r != 3; ++r)
    {
        (void)yepCore_Multiply_V64fV64f_V64f(simplex_batch.row[r], bary1, prod, n_comp);
        (void)yepCore_Add_V64fV64f_V64f(prod, simp1, simp1, n_comp);
    }
    for (r = 0; r != 3; ++r)
    {
        (void)yepCore_Multiply_V64fV64f_V64f(simplex_batch.row[r], bary2, prod, n_comp);
        (void)yepCore_Add_V64fV64f_V64f(prod, simp2, n_comp);
    }
        
    (void)yepCore_Subtract_V64fV64f_V64f(simp1, simp2, diff, n_comp);

    (void)yepCore_Multiply_V64fV64f_V64f(diff, diff, diffsq, n_comp);

    /* Now sum horizontally */
    unsigned p;
    double *pt = diffsq;
    for (p = 0; p != n_points; ++p)
    {
        sqdist[p] = pt[0] + pt[1] + pt[2];
        pt += 3;
    }

    free(prod);
    free(simp1);
    free(simp2);
    free(diff);
    free(diffsq);
}
#endif



double square_dist3(double *a, double *b)
{
    return 
        gsl_pow_2(a[0] - b[0])
        + gsl_pow_2(a[1] - b[1])
        + gsl_pow_2(a[2] - b[2]);
}

double dist3(double *a, double *b)
{
    return sqrt(square_dist3(a, b));
}


double determinant4(double *m)
{
    double d =
        + m[0] * (
                  + m[5] * (m[10] * m[15] - m[14] * m[11])
                  - m[6] * (m[9]  * m[15] - m[11] * m[13])
                  + m[7] * (m[9]  * m[14] - m[10] * m[13])
                  )
        - m[1] * (
                  + m[4] * (m[10] * m[15] - m[14] * m[11])
                  - m[6] * (m[8]  * m[15] - m[11] * m[12])
                  + m[7] * (m[8]  * m[14] - m[10] * m[12])
                  )
        + m[2] * (
                  + m[4] * (m[9] * m[15] - m[11] * m[13])
                  - m[5] * (m[8] * m[15] - m[11] * m[12])
                  + m[7] * (m[8] * m[13] - m[9]  * m[12])
                  )
        - m[3] * (
                  + m[4] * (m[9] * m[14] - m[10] * m[13])
                  - m[5] * (m[8] * m[14] - m[10] * m[12])
                  + m[6] * (m[8] * m[13] - m[9]  * m[12])
                  );
    return d;
}

/* see http://steve.hollasch.net/cgindex/geometry/ptintet.html */
unsigned inside_simplex(double *x)
{
    double dval[4] = { -1, -1, -1, -1 };
    double det[4][4];

    unsigned d;
    for (d = 0; d != 4; ++d)
    {
        memcpy(det, det_template, sizeof(det_template));
        det[d][0] = x[0];
        det[d][1] = x[1];
        det[d][2] = x[2];
        dval[d] = determinant4((double *)det);
        if (dval[d] < 0) break;
    }
    /* sign of det0 is positive (see above) */
    return dval[0] >= 0 && dval[1] >= 0 && dval[2] >= 0 && dval[3] >= 0;
}

