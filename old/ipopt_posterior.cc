// Copyright (C) 2005, 2006 International Business Machines and others.
// All Rights Reserved.
// This code is published under the Common Public License.
//
// $Id: hs071_nlp.cpp 1324 2008-09-16 14:19:26Z andreasw $
//
// Authors:  Carl Laird, Andreas Waechter     IBM    2005-08-16

#include "henry/ipopt_posterior.h"
#include "henry/error_estimate.h"

// for printf
#ifdef HAVE_CSTDIO
# include <cstdio>
#else
# ifdef HAVE_STDIO_H
#  include <stdio.h>
# else
#  error "don't have header file for stdio"
# endif
#endif

using namespace Ipopt;

// constructor
IpoptPosterior::IpoptPosterior(ErrorEstimate const& error_estimate,
                               size_t ndim) :
    ErrorEstimate(error_estimate),
    ndim(ndim)
{}

//destructor
IpoptPosterior::~IpoptPosterior()
{}

// returns the size of the problem
bool IpoptPosterior::get_nlp_info(Index& n, Index& m, Index& nnz_jac_g,
                                  Index& nnz_h_lag, IndexStyleEnum& index_style)
{
    // The problem described in IpoptPosterior.h has 4 variables, x[0] through x[3]
    n = ndim;

    // one equality constraint
    m = 1;

    // in this example the jacobian is dense and contains 8 nonzeros
    nnz_jac_g = ndim;

    // the hessian is also dense and has 16 total nonzeros, but we
    // only need the lower left corner (since it is symmetric)
    nnz_h_lag = ndim;

    // use the C style indexing (0-based)
    index_style = TNLP::C_STYLE;

    return true;
}

// returns the variable bounds
bool IpoptPosterior::get_bounds_info(Index n, Number* x_l, Number* x_u,
                                     Index m, Number* g_l, Number* g_u)
{
    // here, the n and m we gave IPOPT in get_nlp_info are passed back to us.
    // If desired, we could assert to make sure they are what we think they are.
    assert(n == ndim);
    assert(m == 1);

    // the variables have lower bounds of 0
    for (Index i=0; i<ndim; i++) {
        x_l[i] = 0.0;
    }

    // the variables have upper bounds of 1
    for (Index i=0; i<ndim; i++) {
        x_u[i] = 1.0;
    }

    //the first constraint expresses the normalization condition, A+C+G+T = 1
    if (ndim == 3)
    {
        g_l[0] = 0.0;
        g_u[0] = 1.0;
    }
    else
    {
        g_l[0] = 1.0;
        g_u[0] = 1.0;
    }

    return true;
}

// returns the initial point for the problem
bool IpoptPosterior::get_starting_point(Index n, bool init_x, Number* x,
                                        bool init_z, Number* z_L, Number* z_U,
                                        Index m, bool init_lambda,
                                        Number* lambda)
{
    // Here, we assume we only have starting values for x, if you code
    // your own NLP, you can provide starting values for the dual variables
    // if you wish
    assert(init_x == true);
    assert(init_z == false);
    assert(init_lambda == false);

    // initialize to the given starting point
    x[0] = 0.25;
    x[1] = 0.25;
    x[2] = 0.25;

    if (ndim == 4)
    {
        x[3] = 0.25;
    }

    return true;
}

// returns the value of the objective function
bool IpoptPosterior::eval_f(Index n, const Number* x, bool new_x, Number& obj_value)
{
    assert(n == ndim);

    if (ndim == 3)
    {
        Number xa[4];
        std::copy(x, x + 3, xa);
        xa[3] = 1.0 - x[0] - x[1] - x[2];
        obj_value = -1.0 * this->Log2Posterior(xa);
    }
    else
    {
        obj_value = -1.0 * this->Log2Posterior(x);
    }

    return true;
}

// return the gradient of the objective function grad_{x} f(x)
bool IpoptPosterior::eval_grad_f(Index n, const Number* x, bool new_x, Number* grad_f)
{
    assert(n == ndim);

    this->log2_posterior_gradient(x, grad_f);
    for (size_t d = 0; d != ndim; ++d)
    {
        grad_f[d] *= -1.0;
    }

    return true;
}


// return the value of the constraints: g(x)
bool IpoptPosterior::eval_g(Index n, const Number* x, bool new_x, Index m, Number* g)
{
    assert(n == ndim);
    assert(m == 1);

    if (ndim == 4)
    {
        g[0] = x[0] + x[1] + x[2] + x[3];
    }
    else
    {
        g[0] = x[0] + x[1] + x[2];
    }

    return true;
}

// return the structure or values of the jacobian
bool IpoptPosterior::eval_jac_g(Index n, const Number* x, bool new_x,
                                Index m, Index nele_jac, Index* iRow, Index *jCol,
                                Number* values)
{
    if (values == NULL) {
        // return the structure of the jacobian, a single row of four
        // non-zero numbers
        iRow[0] = 0;
        jCol[0] = 0;
        iRow[1] = 0;
        jCol[1] = 1;
        iRow[2] = 0;
        jCol[2] = 2;

        if (ndim == 4)
        {
            iRow[3] = 0;
            jCol[3] = 3;
        }
    }
    else {
        // return the values of the jacobian of the constraints
        // the constraint is a simple sum.  the derivs are constant
        
        values[0] = 1.0;
        values[1] = 1.0;
        values[2] = 1.0;

        if (ndim == 4)
        {
            values[3] = 1.0;
        }

    }

    return true;
}

//return the structure or values of the hessian
bool IpoptPosterior::eval_h(Index n, const Number* x, bool new_x,
                            Number obj_factor, Index m, const Number* lambda,
                            bool new_lambda, Index nele_hess, Index* iRow,
                            Index* jCol, Number* values)
{
    return false;
}


void IpoptPosterior::finalize_solution(SolverReturn status,
                                       Index n, const Number* x, const Number* z_L, const Number* z_U,
                                       Index m, const Number* g, const Number* lambda,
                                       Number obj_value,
                                       const IpoptData* ip_data,
                                       IpoptCalculatedQuantities* ip_cq)
{
    // here is where we would store the solution to variables, or
    // write to a file, etc so we could use the solution.

    // For this example, we write the solution to the console
    printf("\n\nSolution of the primal variables, x\n");
    for (Index i=0; i<n; i++) {
        printf("x[%d] = %e\n", i, x[i]);
    }

    printf("\n\nSolution of the bound multipliers, z_L and z_U\n");
    for (Index i=0; i<n; i++) {
        printf("z_L[%d] = %e\n", i, z_L[i]);
    }
    for (Index i=0; i<n; i++) {
        printf("z_U[%d] = %e\n", i, z_U[i]);
    }

    printf("\n\nObjective value\n");
    printf("f(x*) = %e\n", obj_value);

    printf("\nFinal value of the constraints:\n");
    for (Index i=0; i<m ;i++) {
        printf("g(%d) = %e\n", i, g[i]);
    }
}
