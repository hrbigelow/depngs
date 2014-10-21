#ifndef _METROPOLIS_H
#define _METROPOLIS_H

#include <gsl/gsl_randist.h>

#include "defs.h"

class ErrorEstimate;
class Dirichlet;

struct mh_metric;

class Metropolis
{
 public:

    ErrorEstimate *integrand;
    Dirichlet *proposal;

    gsl_rng *randgen;


    

    double current_point[NUM_NUCS];
    double *sample_points;
    size_t num_points;

    Metropolis(ErrorEstimate *posterior,
               Dirichlet *dirichlet,
               size_t total_sample_points);

    ~Metropolis();

    void propose(double *z_star);

    void sample(size_t const num_samples,
                size_t const burn_in,
                size_t const every_nth,
                double *initial_point,
                double *proposal_mean,
                double *proposal_variance,
                double *alt_sample_points);
    
    double step(struct mh_metric **z_tau, struct mh_metric **z_nxt);

    void tune_proposal(size_t num_points,
                       double max_tuning_iterations,
                       size_t averaging_window_size,
                       size_t front_back_window_distance,
                       size_t autocor_max_offset);
    
};


void print_mean_covariance(FILE *fh, 
                           const double *mean, 
                           const double *covariance, 
                           size_t ndim);


#endif // _METROPOLIS_H
