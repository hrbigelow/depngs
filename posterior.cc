#include "posterior.h"
#include "error_estimate.h"

#include "stats_tools.h"

//currently, initializes mode point and log_scaling_factor
void Posterior::initialize(double mode_tolerance, size_t max_function_evals,
                           double const* initial_point,
                           bool verbose)
{

    this->ee->find_mode_point(1e-20, max_function_evals, 
                              initial_point,
                              this->zero_boundary,
                              verbose,
                              this->mode_point);
    
    this->log_scaling_factor = 
        this->ee->log_likelihood(this->mode_point)
        + this->ee->log_dirichlet_prior(this->mode_point);
}


double Posterior::pdf(double const* x)
{
    double y;
//     static size_t recall_count = 0;
    // if (this->get_last_call(x, &y))
    // {
    //     return y;
    // }

    double xx[4];
    std::copy(x, x+3, xx);
    xx[3] = 1.0 - x[0] - x[1] - x[2];
    if (normalized(xx, 4, 1e-10) && all_positive(xx, 4))
    {
        y = ee->ScaledPosterior(xx, this->log_scaling_factor);
    }
    else
    {
        y = 0.0;
    }
    // this->store_call(x, y);
    return y;
}


double Posterior::log_pdf(double const* x)
{

    double y;

    double xx[4];
    std::copy(x, x+4, xx);
    if (this->ndim == 3)
    {
        xx[3] = std::max(0.0, 1.0 - x[0] - x[1] - x[2]);
    }

    if (normalized(xx, 4, 1e-10) && all_positive(xx, 4))
    {
        y = ee->log_likelihood(xx)
            + ee->log_dirichlet_prior(xx);
    }
    else
    {
        y = -DBL_MAX;
    }
    return y;
}


// double Posterior::log2_pdf(double const*x)
// {
//     double ret = log2_likelihood(ee, x);
//         + ee->log_dirichlet_prior(x);
//     return ret;
// }
