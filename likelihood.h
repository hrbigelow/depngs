#ifndef _LIKELIHOOD_H
#define _LIKELIHOOD_H

/* Computes log2(prod_i(dotp(point_p, cpd_i)^ct_i)) for P points (in
   4D) given the cpd (conditional probability distribution) and counts
   ct for each of the I terms.

   The approach is to batch each step of the calculation across all
   points, rather than perform each step in order for each point.
   This allows vectorization.

   log2(prod_i(dotp(point_p, cpd_i)^ct_i))
   
   = log2(prod_i(q_pi^ct_i)) -- let q_pi = dotp(point_p, cpd_i)
   
   = log2(prod_i((m_pi * 2^e_pi)^ct_i)) -- decompose to mantissa and exponent
   
   = log2(prod_i(m_pi^ct_i)) + log2(prod_i(2^(e_pi * ct_i)))

   = log2(prod_i(m_pi^ct_i)) + sum_i(e_pi * ct_i)

   Notice that the m_pi's are all in [0.5, 1), and so if sum(ct_i) <
   1022, then we can guarantee that the term prod_i(m_pi^ct_i) will
   not underflow.
   
   For each of the I terms:

   1.  Calculate all dot products, store in dotp.
   2.  Extract mantissas and exponents from dotp, store in mnt and exp
   3.  Exponentiate mantissas with ct_i, store in cmnt
   4.  Multiply exponents by ct_i, store in cexp
   5.  Add cexp to scexp
   6.  If (do_scale) zero out exponents of cmnt, adding them to scexp
   7.  Multiply pcmnt by cmnt
   8.  compute log2(pcmnt) * scexp, storing in ll (log2 likelihood)
*/

#include <stddef.h>

struct likelihood_buf {
    double *buf; /* supports memory for all following pointers */
    double *dotp; /* dot products */
    double *mnt, *cmnt, *pcmnt; /* mantissas, mnt^ct, and product of
                                   mnt^ct */
    double *exp, *cexp, *scexp; /* extracted exponents of dotp
                                   numbers. cexp is exp * ct, scexp is
                                   sum of cexp over each of the terms
                                   in the polynomial */
    size_t n_elem; /* number of elements of space in this workspace */
};

void log2_likelihood(struct likelihood_buf bufs,
                     double *points,
                     size_t n_points,
                     double *cpd,
                     unsigned int *counts,
                     size_t n_terms,
                     double *l2l);

#endif /* _LIKELIHOOD_H */
