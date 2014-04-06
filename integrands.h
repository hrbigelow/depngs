#ifndef _INTEGRANDS_H
#define _INTEGRANDS_H

#include "sampling.h"
#include "tools.h"

#include <gsl/gsl_errno.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_blas.h>


class Integrand
{
protected:
    size_t ndim;
    
private:
    double * last_point;
    double last_value;
    
    
public:

    bool may_underflow;

    Integrand(size_t _ndim, bool _may_underflow) : 
        ndim(_ndim), may_underflow(_may_underflow)
    { 
        last_point = new double[ndim];
    }
    ~Integrand()
    {
        delete last_point;
    }
    void store_call(double const* x, double y);
    bool get_last_call(double const* x, double * y) const;

    virtual double pdf(double const* x) = 0;
    virtual double log_pdf(double const* x) = 0;

};

class SamplingFunction : public Integrand
{
public:
    SamplingFunction(size_t _ndim, bool _may_underflow) : Integrand(_ndim, _may_underflow) { }
    virtual void sample(double * x) const = 0;
    virtual void sample_conditioned(double const* x_tau,
                                    double * x_star) = 0;
};




#endif // _INTEGRANDS_H
