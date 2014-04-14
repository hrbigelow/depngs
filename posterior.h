#ifndef _POSTERIOR_H
#define _POSTERIOR_H

#include <cstddef>

class ErrorEstimate;

class Posterior
{

 public:

    ErrorEstimate * ee;
    size_t ndim;
    bool may_underflow;
    
    double log_scaling_factor;
    
    double mode_point[4];
    bool zero_boundary[4];
    
 Posterior(ErrorEstimate * ee, size_t ndim, bool may_underflow) :
    ee(ee), ndim(ndim), may_underflow(may_underflow)
    {
        
    }
    
    void initialize(double mode_tolerance, size_t max_function_evals, 
                    double const* initial_point, bool verbose);

    double pdf(double const* x);
    double log_pdf(double const* x);

};

#endif // _POSTERIOR_H
