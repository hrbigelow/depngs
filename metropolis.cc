#include "metropolis.h"
#include "dirichlet.h"
#include "error_estimate.h"
#include "stats_tools.h"

#include <cassert>
#include <cmath>
#include <algorithm>
#include <ctime>
#include <vector>
#include <string.h>

#include "sampling.h"
#include <gsl/gsl_statistics_double.h>

Metropolis::Metropolis(ErrorEstimate *integrand,
                       Dirichlet *proposal,
                       size_t num_points) 
    : integrand(integrand), 
      proposal(proposal), 
      num_points(num_points)
{ 
    sample_points = new double[NUM_NUCS * num_points];
    this->randgen = gsl_rng_alloc(gsl_rng_taus);
    gsl_rng_set(this->randgen, 0);
}


Metropolis::~Metropolis()
{
    delete sample_points;
    gsl_rng_free(this->randgen);
}


//propose a new point z_star from a starting point z_tau
//z_tau is unused here since Metropolis-Hastings is using one
//unchanging proposal distribution
void Metropolis::propose(double *z_star)
{
    this->proposal->sample(z_star);
}

struct mh_metric
{
    double x[NUM_NUCS];
    double log_ratio_ivp; // integrand vs. proposal
};


double Metropolis::accept(const double *z_star, const double *z_tau)
{
    double ratio;

    double y_tau_log = this->integrand->log_pdf(z_tau);
    double y_star_log = this->integrand->log_pdf(z_star);
    double proposal_y_tau_log = this->proposal->log_pdf(z_tau);
    double proposal_y_star_log = this->proposal->log_pdf(z_star);
    
    assert(! isnan(y_tau_log));
    assert(! isnan(proposal_y_star_log));
    
    double numerator_log = y_star_log - proposal_y_star_log;
    double denominator_log = y_tau_log - proposal_y_tau_log;
    
    //master formula for metropolis hastings
    double ratio_log = numerator_log - denominator_log;
    
    assert(! isnan(ratio_log));
    
    //may be infinite.  this is okay.
    ratio = exp(ratio_log);
    return ratio;

}



void Metropolis::set_current_point(const double *point)
{
    std::copy(point, point + NUM_NUCS, current_point);
}


// pre-condition: **z_tau is the current point, **z_nxt is scratch
// space.  post-condition: **z_tau is the updated current point, which
// is either the same as it was before if the proposal was rejected,
// or it is new.  **z_nxt is again valid scratch space. returns the
// proposal value that was used in the accept/reject logic.
double Metropolis::step(mh_metric **z_tau, mh_metric **z_nxt)
{

    this->propose((*z_nxt)->x);
    (*z_nxt)->log_ratio_ivp = 
        this->integrand->log_pdf((*z_nxt)->x)
        - this->proposal->log_pdf((*z_nxt)->x);

    double ratio_log = 
        (*z_nxt)->log_ratio_ivp 
        - (*z_tau)->log_ratio_ivp;

    assert(! isnan(ratio_log));

    // may be infinite.  this is okay
    double proposal_value = exp(ratio_log);

    mh_metric *tmp;
    if (proposal_value > gsl_rng_uniform(this->randgen))
    {
        tmp = *z_tau;
        *z_tau = *z_nxt;
        *z_nxt = tmp;
    }

    return proposal_value;

}



//do we need to record the y values here?
// here, z_star is scratch space.  
// z_tau is unchanged.
// z_tau_next is assumed to be
// uninitialized, and will be initialized after this call, either
// with the new proposed value (z_star), or the previous value (z_tau)
double Metropolis::step_old(const double *z_tau, double *z_tau_next)
{

    //update the markov chain at each step.
    double proposal_value;

    double z_star[NUM_NUCS];
    double uniform;

    //propose a z_star out of the blue
    this->propose(z_star);

    assert(normalized(z_star, NUM_NUCS, 1e-10));

    //from the integrand values, calculate and acceptance score
    proposal_value = this->accept(z_star, z_tau);

    uniform = gsl_rng_uniform(this->randgen);
            
    if (proposal_value > uniform)
    {
        //new jump is accepted; assign
        std::copy(z_star, z_star + NUM_NUCS, z_tau_next);
    }

    else
    {
        //old position kept; assign
        std::copy(z_tau, z_tau + NUM_NUCS, z_tau_next);
    }

    return proposal_value;

}



// sample from the integrand using the Metropolis algorithm.
// populates the sample buffer.  if alt_sample_points are provided,
// use those instead of internal.  if proposal_mean and
// proposal_variance are non-null, 
void 
Metropolis::sample(size_t num_samples_to_take,
                   size_t burn_in,
                   size_t every_nth,
                   double *proposal_mean,
                   double *proposal_variance,
                   double *alt_sample_points)
{
    if (num_samples_to_take > this->num_points)
    {
        fprintf(stderr, "Metropolis::sample:  %Zu samples requested exceeds "
                "%Zu samples allocated for this Metropolis instance"
                "Please call the constructor with larger number of samples", 
                num_samples_to_take, this->num_points);
        exit(1);
    }

    double 
        *point = alt_sample_points ? alt_sample_points : this->sample_points,
        *point_end = point + (num_samples_to_take * NUM_NUCS);

    //update the markov chain at each step.
    srand(time(NULL));
    double proposal_value;
    size_t step_count = 0;

    // scratch space
    struct mh_metric m1, m2;

    mh_metric *z_tau = &m1, *z_star = &m2;

    // std::copy(this->current_point, this->current_point + NUM_NUCS, z_tau->x);
    memcpy(z_tau->x, this->current_point, sizeof(double) * NUM_NUCS);
    z_tau->log_ratio_ivp = 
        this->integrand->log_pdf(z_tau->x)
        - this->proposal->log_pdf(z_tau->x);


    double *proposal_values = NULL;
    if (proposal_mean)
        proposal_values = new double[num_samples_to_take];

    while (point != point_end)
    {
        proposal_value = this->step(&z_tau, &z_star);

        if (proposal_mean && step_count < num_samples_to_take)
            proposal_values[step_count] = 
                std::min(1.0, isinf(proposal_value) ? DBL_MAX : proposal_value);

        if (step_count > burn_in && step_count % every_nth == 0)
        {
            memcpy(point, z_tau->x, sizeof(double) * NUM_NUCS);
            // std::copy(z_tau->x, z_tau->x + NUM_NUCS, point);
            point += NUM_NUCS;
        }

        ++step_count;

    }

    if (proposal_values)
    {
        *proposal_mean = gsl_stats_mean(proposal_values, 1, num_samples_to_take);
        *proposal_variance = gsl_stats_variance_m(proposal_values, 1, num_samples_to_take, *proposal_mean);
    }

    memcpy(this->current_point, z_tau->x, sizeof(double) * NUM_NUCS);
    // this->set_current_point(z_tau);

    if (proposal_values) 
        delete proposal_values;
}





//  other support functions
void print_mean_covariance(FILE *fh, const double *mean, 
                           const double *covariance)
{
    for (size_t d1 = 0; d1 != NUM_NUCS; ++d1)
    {
        fprintf(fh, "%20.18g\t", mean[d1]);

        for (size_t d2 = 0; d2 != NUM_NUCS; ++d2)
            fprintf(fh, "\t%20.18g", covariance[d1 * NUM_NUCS + d2]);

        fprintf(fh, "\n");
    }
}
