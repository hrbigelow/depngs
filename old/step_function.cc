#include "henry/step_function.h"


StepFunction::StepFunction(double const* _low_bounds,
                           double const* _hi_bounds,
                           size_t ndim, bool _is_log2_integrand) 
    : AnalyticalIntegrand(ndim, _is_log2_integrand)
{
    this->low_bounds = new double[this->ndim];
    this->hi_bounds = new double[this->ndim];
    std::copy(_low_bounds, _low_bounds + this->ndim, this->low_bounds);
    std::copy(_hi_bounds, _hi_bounds + this->ndim, this->hi_bounds);
}


StepFunction::~StepFunction()
{
    delete this->low_bounds;
    delete this->hi_bounds;
}


double StepFunction::pdf(double const* x)
{
    bool in_step = true;
    for (size_t d = 0; d != this->ndim; ++d)
    {
        in_step = in_step 
            && x[d] >= this->low_bounds[d] 
            && x[d] <= this->hi_bounds[d];
    }
    return in_step ? 1.0 : 0.0;
}


double StepFunction::marginal_cdf(double const x, size_t const marg_dim) const
{
    double l = this->low_bounds[marg_dim];
    double h = this->hi_bounds[marg_dim];
    double w = h - l;
    if (x < l)
    {
        return 0.0;
    }
    else if (x < h)
    {
        return (x - l) / w;
    }
    else
    {
        return 1.0;
    }
}


double StepFunction::inv_marginal_cdf(double const p, size_t const marg_dim) const
{
    double l = this->low_bounds[marg_dim];
    double h = this->hi_bounds[marg_dim];
    double w = h - l;
    if (p == 0)
    {
        return l;
    }
    else if (p < 1.0)
    {
        return l + w * p;
    }
    else
    {
        return h;
    }
}
