#ifndef _DIRICHLET_H
#define _DIRICHLET_H

#include <cstddef>
#include <gsl/gsl_rng.h>

#include "defs.h"

class Dirichlet
{
public:

    gsl_rng *seed;
    double alpha[NUM_NUCS], alpha0;
    double Zlog2;

    Dirichlet();
    ~Dirichlet();

    void update(const double *_alpha);

    void set_alphas_from_mode(const double *mode);
    void set_alphas_from_mean(const double *mean);

    void set_alphas_from_mode_or_bound(const double *mode,
                                       const double *alpha_lower_bound,
                                       bool const* is_zero_boundary);

    void set_alphas_from_mean_or_bound(const double *mean,
                                       const double *lower_bound);

    void lower_bound_alphas(const double *lower_bound);

    double pdf(const double *x);
    double log_pdf(const double *x);

    double log_pdf_4d(const double *x);
    void sample(double * x) const;
};


double
ran_dirichlet_lnpdf(const size_t K,
                    const double alpha[], const double theta[]);

#endif // _DIRICHLET_H
