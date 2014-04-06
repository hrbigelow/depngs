#include "integrands.h"

#include <gsl/gsl_cdf.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_permutation.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_sf_log.h>

#include "stats_tools.h"
#include "error_estimate.h"

// not thread safe
void Integrand::store_call(double const* x, double y)
{
    std::copy(x, x + this->ndim, this->last_point);
    this->last_value = y;
}    


bool Integrand::get_last_call(double const* x, double * y) const
{
    *y = this->last_value;
    return std::equal(this->last_point, this->last_point + this->ndim, x);
}

void ignore_underflow_handler(char const* /* reason */,
                              char const* /* file */,
                              int /* line */,
                              int /* gsl_errno */) { }
