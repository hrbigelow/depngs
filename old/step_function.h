#ifndef _STEP_FUNCTION_H
#define _STEP_FUNCTION_H

#include "henry/integrands.h"

class StepFunction : public AnalyticalIntegrand
{
    double * low_bounds;
    double * hi_bounds;

public:

    StepFunction(double const* low_bounds,
                 double const* hi_bounds,
                 size_t ndim,
                 bool _is_log2_integrand);

    ~StepFunction();

    REAL pdf(double const* x);
    void sample(double * x) const;
    double marginal_cdf(double const xcoord, size_t const marg_dim) const;
    double inv_marginal_cdf(double const p, size_t const marg_dim) const;
};


#endif // _STEP_FUNCTION_H
