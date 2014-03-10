#include "error_estimate.h"

#include <cstdio>
#include <cstdlib>
#include <cassert>
#include <cmath>
#include <cstring>
#include <limits>
#include <utility>
#include <map>
#include <algorithm>
#include <set>
#include <vector>
#include <cctype>

#include <gsl/gsl_vector.h>
#include <gsl/gsl_multimin.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_sf_log.h>
#include <gsl/gsl_sf_exp.h>
#include <gmp.h>

#include "hilbert.h"
#include "dirichlet.h"

#include "tools.h"
#include "sampling.h"
#include "stats_tools.h"
#include "simulation.h"
#include "transformation.h"
#include "nucleotide_stats.h"


/*
  In the future, we will be extending the Posterior analysis to the
more general case of the data being not a (base, qual) pair, but a set
of four estimated probabilities of the underlying 'picked' base p,
P(p).  These will be based on the image evidence which we summarize as
P(p|i).  
 */


const int ErrorEstimate::NBASES;

REAL expansion_rows_global[3][3] = { 
    { 1., 1., 1. }, 
    { 0., 2., 1. }, 
    { 0., 0., 3. } 
};

REAL const expansion_determinant = 6.0;
REAL const contraction_determinint = 1.0 / 6.0;

REAL contraction_rows_global[3][3] = { 
    { 1., -0.5, -1.0/6.0 }, 
    { 0., 0.5, -1.0/6.0 }, 
    { 0., 0., 1.0/3.0 } 
};


ErrorEstimate::ErrorEstimate()
{
    //create matrices for transformations
    std::fill(this->composition_prior_alphas,
              this->composition_prior_alphas + 4, 1.0);
    this->uniform_prior = true;
    this->num_discrete_priors = 0;
    this->log_discrete_prior_dist = NULL;
    this->prior_type = ErrorEstimate::CONTINUOUS;

    memcpy(expansion_rows, expansion_rows_global, 9 * sizeof(REAL));
    memcpy(contraction_rows, contraction_rows_global, 9 * sizeof(REAL));
}


void ErrorEstimate::set_composition_prior_alphas(double const* alphas)
{
    std::copy(alphas, alphas + 4, this->composition_prior_alphas);
    this->uniform_prior = 
        alphas[0] == 1.0
        && alphas[1] == 1.0
        && alphas[2] == 1.0
        && alphas[3] == 1.0;
}


void ErrorEstimate::set_discrete_prior_dist(double const* prior_points_flat,
                                            double const* prior_dist, 
                                            size_t num_priors)
{
    this->num_discrete_priors = num_priors;
    this->log_discrete_prior_dist = new double[this->num_discrete_priors];
    std::copy(prior_dist, prior_dist + num_priors, this->log_discrete_prior_dist);

    double const* prior_point;
    for (size_t i = 0; i != num_priors; ++i)
    {
        prior_point = prior_points_flat + (i * 4);
        this->discrete_prior_index[prior_point] = i;

        this->log_discrete_prior_dist[i] = 
            log(this->log_discrete_prior_dist[i]);
    }
    this->prior_type = ErrorEstimate::DISCRETE;
}


ErrorEstimate::~ErrorEstimate()
{
    if (this->log_discrete_prior_dist != NULL)
    {
        delete this->log_discrete_prior_dist;
    }
}


//calculates P(Obs,sample_comp) as sum_fb { P(sample_comp)P(fb|sample_comp)P(Obs|fb) }
//can never be zero
// !!! Note:  This can be optimized as a vectorized dot product
/*
double ErrorEstimate::single_observation(double const* sample_comp,
                                         size_t di) const
{
    assert(di < this->model_params->num_distinct_data);

    double return_val = 
        sample_comp[0] * this->model_params->founder_base_likelihood[0][di]
        + sample_comp[1] * this->model_params->founder_base_likelihood[1][di]
        + sample_comp[2] * this->model_params->founder_base_likelihood[2][di]
        + sample_comp[3] * this->model_params->founder_base_likelihood[3][di];
    
    if (isnan(return_val))
    {
        assert(! isnan(return_val));
        assert(! isinf(return_val));
    }
    //assert(return_val > 0.0); //could be 0 if we have a case of infinite quality
    return return_val;

}
*/


// calculates d/dC P(I_b,C) as sum_b { P(C)P(b|C)P(I|b) }
// What?!  This doesn't seem to depend on the actual sample composition
// 
 /*
double ErrorEstimate::single_observation_gradient(double const* sample_composition,
                                                  size_t datum_index,
                                                  size_t deriv_dimension) const
{
    // assert(datum_index < this->model_params->num_distinct_data);
    return this->model_params->founder_base_likelihood[deriv_dimension][datum_index];
}
 */


//calculate d/dC ( log P(C,I_1,...,I_D) )

/*
void ErrorEstimate::log_likelihood_gradient(double const* sample_composition,
                                            double * gradient) const
{

    size_t ndim = 4;

    std::fill(gradient, gradient + ndim, 0.0);
    for (size_t d = 0; d != ndim; ++d)
    {
        //sum_g(frac{1}{ln(2)P(I|C)} P(b|C))
        for (size_t raw_idx = 0; 
             raw_idx != this->locus_data->num_distinct_data; ++raw_idx)
        {
            size_t datum_index = this->locus_data->stats_index[raw_idx];
            double count = this->locus_data->raw_counts[raw_idx];
            if (count == 0)
            {
                continue;
            }
            gradient[d] +=
                (count * this->single_observation_gradient(sample_composition, datum_index, d))
                / this->single_observation(sample_composition, datum_index);
        }
//         gradient[d] /= M_LOG2E;
//         gradient[d] += log_comp_prior_gradient(sample_composition, d);
    }
}
*/

//calculate d/dC ( log P(C,I_1,...,I_D) )
void ErrorEstimate::log_likelihood_gradient(double const* comp,
                                            double * gradient) const
{

    std::fill(gradient, gradient + 4, 0.0);

    double * l = this->locus_data->fbqs_cpd;
    double * l_end = l + (this->locus_data->num_data * 4);
    unsigned long * lc = this->locus_data->raw_counts;

    //sum_g(frac{1}{ln(2)P(I|C)} P(b|C))

    // iterate over each bqs category
    for (; l != l_end; l += 4, lc++)
    {
        double so =
            (*l) * comp[0] + (*(l+1)) * comp[1] + (*(l+2)) * comp[2] + (*(l+3)) * comp[3];
    
        gradient[0] += (*lc) * (*l) / so;
        gradient[1] += (*lc) * (*(l+1)) / so;
        gradient[2] += (*lc) * (*(l+2)) / so;
        gradient[3] += (*lc) * (*(l+3)) / so;
    }
}



double ErrorEstimate::log_discrete_prior(size_t sample_point_index) const
{
    return this->log_discrete_prior_dist[sample_point_index];
}


double ErrorEstimate::log_dirichlet_prior(double const* sample_composition) const
{
    double retval;
    if (this->uniform_prior)
    {
        retval = 0.0;
    }
    else
    {
        retval = 
            Transformation::log_dirichlet(this->composition_prior_alphas,
                                          sample_composition);
    }
    return retval;
    //return (isnan(retval) || isinf(retval)) ? FLT_MAX : retval;
}


/*
  Derivation:
  d/dC_i[log(~Dir(C))] 
  = d/dC_i[log(e) ln(~Dir(C))]
  = log(e) 1/~Dir(C) d/dC_i[~Dir(C)]
  = log(e) (alpha_i - 1) / C_i   // because of cancellation of terms
 */
// double ErrorEstimate::log_composition_prior_gradient(double const* sample_composition,
//                                                       size_t deriv_dimension) const
// {
//     double retval =
//         M_LOG2E
//         * (this->composition_prior_alphas[deriv_dimension] - 1.0)
//         / sample_composition[deriv_dimension];

//     return (isnan(retval) || isinf(retval)) ? FLT_MAX : retval;
// }

/*
REAL ErrorEstimate::log_likelihood(double const* sample_composition) const
{
    
    if (! (normalized(sample_composition, 4, 1e-10) &&
           all_positive(sample_composition, 4)))
    {
        fprintf(stderr, "log_likelihood: invalid input.\n");
        exit(2);
    }
                
    REAL sum_log_factors = 0.0;

    for (size_t raw_index = 0; 
         raw_index != this->locus_data->num_distinct_data; ++raw_index)
    {

        double count = this->locus_data->raw_counts[raw_index];

        if (count == 0)
        {
            //no counts of this data, thus no contribution to log_likelihood
            continue;
        }

        size_t datum_index = this->locus_data->stats_index[raw_index];

        sum_log_factors += 
            gsl_sf_log(single_observation(sample_composition, datum_index))
            * static_cast<REAL>(count);
    }
    return sum_log_factors;
}
*/

 // In this new formulation, we use the subset of the model that is packed into
 // the locus itself.  Also, we inline the 'single_observation' function
double ErrorEstimate::log_likelihood(double const* comp) const
{
    
    // !!! This can be taken out at some point
    // if (! (normalized(comp, 4, 1e-10) && all_positive(comp, 4)))
    // {
    //     fprintf(stderr, "log_likelihood: invalid input.\n");
    //     exit(2);
    // }
                
    double sum_log_factors = 0.0;

    double * l = this->locus_data->fbqs_cpd;
    double * l_end = l + (this->locus_data->num_data * 4);
    unsigned long * lc = this->locus_data->raw_counts;

    // iterate over each BQS category
    for (; l != l_end; l += 4, lc++)
    {
        sum_log_factors += (*lc)
            * gsl_sf_log((*l) * comp[0]
                         + (*(l+1)) * comp[1]
                         + (*(l+2)) * comp[2]
                         + (*(l+3)) * comp[3]);
    }   

    // mpf_t term, pterm, prod;
    // mpf_init(term);
    // mpf_init(pterm);
    // mpf_init_set_ui(prod, 1);

    // for (; l != l_end; l += 4, lc++)
    // {
    //     mpf_set_d(term, 
    //               (*l) * comp[0]
    //               + (*(l+1)) * comp[1]
    //               + (*(l+2)) * comp[2]
    //               + (*(l+3)) * comp[3]);

    //     mpf_pow_ui(pterm, term, (*lc));
    //     mpf_mul(prod, prod, pterm);
    // }   
    // long int exp;
    // double mant = mpf_get_d_2exp(&exp, prod);
    // sum_log_factors = gsl_sf_log(mant) + exp * M_LN2;
    // mpf_clear(term);
    // mpf_clear(pterm);
    // mpf_clear(prod);

    // if (sum_log_factors != 0)
    // {
    //     assert(sum_log_factors / sum_log_factors_mpf < 1.000000001);
    //     assert(sum_log_factors_mpf / sum_log_factors < 1.000000001);
    // }

    return sum_log_factors;
}




//     double comp_prior;
//     switch (this->prior_type)
//     {
//     case ErrorEstimate::DISCRETE:
//         {
//             std::map<double const*, size_t>::const_iterator dp_iter = 
//                 this->discrete_prior_index.find(sample_composition);
//             if (dp_iter == this->discrete_prior_index.end())
//             {
//                 fprintf(stderr, "ErrorEstimate::Log2Posterior: using DISCRETE prior"
//                         " but called with unrecognized composition\n");
//                 exit(1);
//             }
//             size_t sample_index = (*dp_iter).second;
//             comp_prior = this->log2_discrete_prior(sample_index);
//         }
//         break;
//     case ErrorEstimate::CONTINUOUS:
//         comp_prior = this->log2_composition_prior(sample_composition);
//         break;
//     }

//     //could be -infinity
//     if (isinf(sum_log_factors))
//     {
//         return -DBL_MAX;
//     }
//     //assert(finite(sum_log_factors) != 0);
//     else
//     {
//         return sum_log_factors + comp_prior;
//     }
// }



//transforms x (bounds x1[0,1], x2[0,1-x1], x3[0,1-x1-x2]) to
//expanded (bounds in unit hypercube).
//points x containing
typedef std::pair<double, size_t> Key;

bool sort_first_desc(Key a, Key b)
{
    return a.first > b.first;
}


void auxiliary_transform(REAL const matrix[3][3], double const x[3], double * transformed)
{
    
    Key key[] = { Key(x[0], 0), Key(x[1], 1), Key(x[2], 2) };
    std::sort(key, key + 3, sort_first_desc);
    
    size_t rev_key[3];
    for (size_t k = 0; k != 3; ++k)
    {
        rev_key[key[k].second] = k;
    }

    //create the expanded coordinate as the matrix * vector product
    //the matrix is the permutation of the original matrix rows

    for (size_t r = 0; r != 3; ++r)
    {
        double sum = 0.0;
        for (size_t c = 0; c != 3; ++c)
        {
            sum += matrix[rev_key[r]][rev_key[c]] * x[c];
        }
        transformed[r] = sum;
    }
}


bool within_hypercube(double const x[3])
{
    return 
        x[0] >= 0.0 && x[0] <= 1.0 &&
        x[1] >= 0.0 && x[1] <= 1.0 &&
        x[2] >= 0.0 && x[2] <= 1.0;
}


bool within_pyramid(double const x[3])
{
    return
        x[0] >= 0.0 &&
        x[1] >= 0.0 &&
        x[2] >= 0.0 &&
        x[0] + x[1] + x[2] <= 1.0;
}


void ErrorEstimate::expand_to_hypercube(double const x[3], double * expanded) const
{
    assert(within_pyramid(x));

    auxiliary_transform(this->expansion_rows, x, expanded);
}

void ErrorEstimate::contract_from_hypercube(double const x[3], double * contracted) const
{
    assert(within_hypercube(x));

    auxiliary_transform(this->contraction_rows, x, contracted);
}


REAL ErrorEstimate::ScaledPosterior(double const* sample_composition,
                                    REAL log_scaling_factor) const
{

    REAL log_scaled =
        this->log_likelihood(sample_composition) 
        + this->log_dirichlet_prior(sample_composition) 
        - log_scaling_factor;
    
    double value = (log_scaled > -1023.0) ? gsl_sf_exp(log_scaled) : 0.0;

    assert(! isnan(value));
    assert(! isinf(value));
    return value;
    
}




//finds mode point of this posterior, returning true on success
//returns number of iterations
size_t ErrorEstimate::find_mode_point(double min_step_size, 
                                      size_t max_iterations,
                                      double const* initial_point,
                                      bool * on_zero_boundary,
                                      bool verbose,
                                      double * mode_point) const
{

    //the minimization is performed in 3 unconstrained dimensions
    size_t sphere_ndim = 3;

    gsl_multimin_function_fdf fdf_minimizer_function;

    Transformation::PassingParams passing_params = {this, 0.0};
    
    fdf_minimizer_function.n = sphere_ndim;
    fdf_minimizer_function.params = & passing_params;
    fdf_minimizer_function.f = &Transformation::log_neg_posterior_value;
    fdf_minimizer_function.df = &Transformation::log_neg_posterior_gradient;
    fdf_minimizer_function.fdf = &Transformation::log_neg_posterior_value_and_gradient;

    gsl_vector *x;
    gsl_multimin_fdfminimizer * first_minimizer =
        gsl_multimin_fdfminimizer_alloc(gsl_multimin_fdfminimizer_vector_bfgs,
                                        sphere_ndim);

    //shirk the initial point.  Just start in one of the four
    //quadrants that gets the highest value.
    // double majority_comp = 0.9;
    // double minority_comp = 0.0333333;
    // double test_comp[4][4];
    double test_values[4];
    // double test_r3_point[sphere_ndim];
    // double best_initial_point[4];

    // precomputed sigmoid values representing 90%, 3.33%, 3.33%, 3.33% base composition.
    double test_r3_point[4][3] = {
        {2.6390567939010965, -0, 3.2958378660048289},
        {2.6390567939010965, -0, -3.2958378660048289},
        {-2.6390584010443274, 3.2958378660048289, -0},
        {-2.6390584010443274, -3.2958378660048289, -0}
    };
    
    for (size_t d = 0; d != 4; ++d)
    {
        // std::fill(test_comp[d], test_comp[d] + 4, minority_comp);
        // test_comp[d][d] = majority_comp;
        // Transformation::composition_to_r3_sigmoid(test_comp[d], 
        //                                           test_r3_point);

        gsl_vector *test_r3_vec = gsl_vector_alloc(sphere_ndim);
        gsl_vector_set(test_r3_vec, 0, test_r3_point[d][0]);
        gsl_vector_set(test_r3_vec, 1, test_r3_point[d][1]);
        gsl_vector_set(test_r3_vec, 2, test_r3_point[d][2]);

        test_values[d] = 
            Transformation::log_neg_posterior_value(test_r3_vec, static_cast<void *>(&passing_params));

        gsl_vector_free(test_r3_vec);
    }
    size_t best_ind = 
        std::distance(test_values, std::min_element(test_values, test_values + 4));

    // std::copy(test_comp[min_test_value_ind], 
    //           test_comp[min_test_value_ind] + 4, 
    //           best_initial_point);

    // double initial_r3_point[sphere_ndim];
    // Transformation::composition_to_r3_sigmoid(best_initial_point, initial_r3_point);

    x = gsl_vector_alloc(sphere_ndim);
    gsl_vector_set(x, 0, test_r3_point[best_ind][0]);
    gsl_vector_set(x, 1, test_r3_point[best_ind][1]);
    gsl_vector_set(x, 2, test_r3_point[best_ind][2]);


    //gsl_vector * numerical_gradient = gsl_vector_alloc(sphere_ndim);
    gsl_vector * last_point = gsl_vector_alloc(sphere_ndim);
    gsl_vector * point_delta = gsl_vector_alloc(sphere_ndim);

    gsl_vector_set_all(last_point, 0.0);

    gsl_multimin_fdfminimizer_set(first_minimizer, &fdf_minimizer_function, x, 0.1, 0.1);

    //double mode_gradient[4][3];

    //double const sigmoid_epsilon = 1e-10;

    Transformation::SigmoidVals sigmoid_vals[3];

    size_t multimin_iteration;
    for (multimin_iteration = 0; multimin_iteration != max_iterations; 
         ++multimin_iteration)
    {

        //copy current point to last point
        gsl_vector_memcpy(last_point, first_minimizer->x);

        //take a step
        gsl_multimin_fdfminimizer_iterate(first_minimizer);

        gsl_vector * gradient = gsl_multimin_fdfminimizer_gradient(first_minimizer);

        double x[3];
        x[0] = gsl_vector_get(first_minimizer->x, 0);
        x[1] = gsl_vector_get(first_minimizer->x, 1);
        x[2] = gsl_vector_get(first_minimizer->x, 2);
        
        Transformation::sigmoid_value_and_gradient(x, sigmoid_vals);
        
        if (verbose)
        {
            //print current stats
            
            Transformation::sigmoid_composition(sigmoid_vals, mode_point);
            printf(
                   "gslfdf1: iter: %Zu, log_neg: %10.8f\t" 
                   "gradient: (%10.8f, %10.8f, %10.8f)\t"
                   "cur_mode: (%10.8f, %10.8f, %10.8f, %10.8f)\t"
                   "r3_point: (%10.8f, %10.8f, %10.8f)\n"
                   ,
                   multimin_iteration,
                   Transformation::log_neg_posterior_value(first_minimizer->x, & passing_params),
                   gsl_vector_get(gradient, 0),
                   gsl_vector_get(gradient, 1),
                   gsl_vector_get(gradient, 2),
                   mode_point[0], mode_point[1], mode_point[2], mode_point[3],
                   gsl_vector_get(first_minimizer->x, 0),
                   gsl_vector_get(first_minimizer->x, 1),
                   gsl_vector_get(first_minimizer->x, 2)
                   );
        }
        
        //check if step is small enough
        gsl_vector_memcpy(point_delta, last_point);
        gsl_vector_sub(point_delta, first_minimizer->x);
        double last_step_size = gsl_blas_dnrm2(point_delta);

        if (gsl_multimin_test_gradient(gradient, 1e-10) == GSL_SUCCESS ||
            (multimin_iteration > 20
             && last_step_size < min_step_size)
            )
        {
            break;
        }

    }

    Transformation::sigmoid_composition(sigmoid_vals, mode_point);
    Transformation::boundary_point(sigmoid_vals, on_zero_boundary);

    gsl_vector_free(last_point);
    gsl_vector_free(point_delta);
    return multimin_iteration;
}
