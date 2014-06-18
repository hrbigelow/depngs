#include "posterior.h"
#include "error_estimate.h"
#include "stats_tools.h"

#include <float.h>

//currently, initializes mode point and log_scaling_factor
void Posterior::initialize(double gradient_tolerance, 
                           size_t max_function_evals,
                           double const* initial_point,
                           bool verbose)
{

    this->ee->find_mode_point(gradient_tolerance, max_function_evals, 
                              initial_point,
                              this->zero_boundary,
                              verbose,
                              this->mode_point);
    
    // this->log_scaling_factor = 
    //     this->ee->log_likelihood(this->mode_point)
    //     + this->ee->log_dirichlet_prior(this->mode_point);
}

/*
  double Posterior::pdf(double const* x)
  {
  double y;
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
  return y;
  }
*/

// Since 'ErrorEstimate' assumes the original point is in 4 dimensions
// and is normalized.  This function attempts to wrap it and allow
// either 3 dimensions or 4.
// But, Metropolis assumes 4 dimensions, while Slice sampling assumes 3.
// As it turns out, slice_sampling.cc just 
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
