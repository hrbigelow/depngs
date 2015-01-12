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


// pre-condition: **z_tau is the current point, **z_nxt is scratch
// space.  post-condition: **z_tau is the updated current point, which
// is either the same as it was before if the proposal was rejected,
// or it is new.  **z_nxt is again valid scratch space. returns the
// proposal value that was used in the accept/reject logic.
/*
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
*/


// compute Log(Post(x) / Prop(x))
double compute_log_ratio_ivp(double *proposal_alpha, double *prior_alpha, 
                             double log_likelihood, double *x)
{
    // Dir(a+b) := Dir(a)Dir(b+1)
    // Post(x) := Dir(a;x)Mult(x)
    // Prop(a+b;x) := Dir(a+b;x) := Dir(a;x)Dir(b+1;x)
    // Post(x) / Prop(x) = Mult(x) / Dir(b+1;x)
    // where a+b = pra (proposal alphas), and
    // a = ipa (integrand prior alphas)
    // So, b+1 = pra - ipa + 1.  Call these 'alpha_residual'
    double alpha_residual[] = { 
        proposal_alpha[0] - prior_alpha[0] + 1, 
        proposal_alpha[1] - prior_alpha[1] + 1,
        proposal_alpha[2] - prior_alpha[2] + 1,
        proposal_alpha[3] - prior_alpha[3] + 1 
    };

    double lr = log_likelihood - gsl_ran_dirichlet_lnpdf(NUM_NUCS, alpha_residual, x);

    assert(! isnan(lr));
    return lr;
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
        compute_log_ratio_ivp(this->proposal->alpha,
                              this->integrand->prior_alpha,
                              this->integrand->log_likelihood((*z_nxt)->x),
                              (*z_nxt)->x);

    assert(! isnan((*z_nxt)->log_ratio_ivp));

    double ratio_log = 
        (*z_nxt)->log_ratio_ivp 
        - (*z_tau)->log_ratio_ivp;

    assert(! isnan(ratio_log));

    // may be infinite.  this is okay
    double proposal_value = exp(ratio_log);

    mh_metric *tmp;

    /* accepting means swapping the values of tau and nxt.  rejecting
       means leaving them the same.  if they are left the same, then */
    if (proposal_value > gsl_rng_uniform(this->randgen))
    {
        tmp = *z_tau;
        *z_tau = *z_nxt;
        *z_nxt = tmp;
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
                   double *initial_point,
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

    memcpy(this->current_point, initial_point, sizeof(double) * NUM_NUCS);

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

    memcpy(z_tau->x, this->current_point, sizeof(double) * NUM_NUCS);

    z_tau->log_ratio_ivp = 
        compute_log_ratio_ivp(this->proposal->alpha,
                              this->integrand->prior_alpha,
                              this->integrand->log_likelihood(z_tau->x),
                              z_tau->x);

    double *proposal_values = NULL;
    if (proposal_mean)
        proposal_values = new double[num_samples_to_take];

    while (point != point_end)
    {
        proposal_value = this->step(&z_tau, &z_star);

        // store the current point as a sample point.  In the previous
        // implmentation, this was erroneously recording the z_tau
        // point before it had been rejected or accepted.  Here, the
        // stored point is the current point after acceptance or rejection.
        if (step_count > burn_in && step_count % every_nth == 0)
        {
            memcpy(point, z_tau->x, sizeof(double) * NUM_NUCS);
            point += NUM_NUCS;
        }

        if (proposal_mean && step_count < num_samples_to_take)
            proposal_values[step_count] = 
                std::min(1.0, isinf(proposal_value) ? DBL_MAX : proposal_value);

        ++step_count;

    }

    if (proposal_values)
    {
        *proposal_mean = gsl_stats_mean(proposal_values, 1, num_samples_to_take);
        *proposal_variance = gsl_stats_variance_m(proposal_values, 1, num_samples_to_take, *proposal_mean);
    }

    memcpy(this->current_point, z_tau->x, sizeof(double) * NUM_NUCS);

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
