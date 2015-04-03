#include "error_estimate.h"

#include <cstring>
#include <cassert>
#include <algorithm>

#include <gsl/gsl_multimin.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_sf_exp.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_randist.h>

#include "transformation.h"
#include "nucleotide_stats.h"
#include "stats_tools.h"

#include <immintrin.h>

ErrorEstimate::ErrorEstimate()
{
    std::fill(this->prior_alpha, this->prior_alpha + 4, 1.0);
}


void ErrorEstimate::set_prior_alphas(const double *alphas)
{
    std::copy(alphas, alphas + 4, this->prior_alpha);
    this->uniform_prior = 
        alphas[0] == 1.0
        && alphas[1] == 1.0
        && alphas[2] == 1.0
        && alphas[3] == 1.0;
}

ErrorEstimate::~ErrorEstimate()
{
}


void ErrorEstimate::log_likelihood_gradient(const double *comp,
                                            double *gradient) const
{
    std::fill(gradient, gradient + 4, 0.0);

    double *l = this->locus_data->fbqs_cpd;
    double *l_end = l + (this->locus_data->num_data * 4);
    unsigned long *lc = this->locus_data->raw_counts;

    //sum_g(frac{1}{ln(2)P(I|C)} P(b|C))

    // iterate over each bqs category
    // __m256d vc, vl, vg, vm;
    
    for (; l != l_end; l += 4, lc++)
    {
        // vc = _mm256_loadu_pd(comp);
        // vl = _mm256_loadu_pd(l);
        // vm = __builtin_ia32_mulpd256(vc, vl);
       
        double so =
            l[0] * comp[0] 
            + l[1] * comp[1] 
            + l[2] * comp[2] 
            + l[3] * comp[3];
    
        gradient[0] += (*lc) * l[0] / so;
        gradient[1] += (*lc) * l[1] / so;
        gradient[2] += (*lc) * l[2] / so;
        gradient[3] += (*lc) * l[3] / so;
    }
}


// In this new formulation, we use the subset of the model that is packed into
// the locus itself.  Also, we inline the 'single_observation' function
/*
double ErrorEstimate::log_likelihood(const double *comp) const
{
    
    double *l = this->locus_data->fbqs_cpd;
    double *l_end = l + (this->locus_data->num_data * 4);
    unsigned long *lc = this->locus_data->raw_counts;
    double sum_log_factors = 0.0;

    // this is the original calculation.  below is an optimized version
    for (; l != l_end; l += 4, lc++)
    {
        double q = (*l) * comp[0] + (*(l+1)) * comp[1] + (*(l+2)) * comp[2] + (*(l+3)) * comp[3];
        sum_log_factors += (*lc) * log2(q);
    }
    
    sum_log_factors *= M_LN2;
    return sum_log_factors;
}
*/



/*
  1. q[i] = dot_prod_i for all q
  2. m[i], e[i] for all q[i] (frexp)
  3. em[i] = m[i] ** ct[i] for all i (these should not overflow or underflow)
  4. ee[i] = e[i] * ct[i] for all i
  5. max_ee = max over i of ee[i]

  In a second loop:
  6. emn[i] = em[i] * 2**(max_ee - ee[i])  (some of these will underflow to zero, silently)
  7. sum_m = sum over i of emn[i]
  8. lsm = log2(sum_m) + max_ee

Finally, transform to base e:

ln(X) = log2(X) / log2(e)


*/


double ErrorEstimate::log_likelihood(const double *x) const
{
    double *cpd_beg = this->locus_data->fbqs_cpd;
    double *cpd_end = cpd_beg + (this->locus_data->num_data * 4);
    unsigned long *ct_beg = this->locus_data->raw_counts;

    double *cpd;
    unsigned long *ct;

    // this is the original calculation.  below is an optimized version
    int s = 0;
    double p = 1.0;
    double lp = 0.0;
    int cti_floor = 0;
    int min_exp = DBL_MIN_EXP;
    for (cpd = cpd_beg, ct = ct_beg; cpd != cpd_end; cpd += 4, ++ct)
    {
        int e, cti = (int)*ct;
        double q = (cpd[0] * x[0]) + (cpd[1] * x[1]) + (cpd[2] * x[2]) + (cpd[3] * x[3]);
        if (cti >= 50 || cti_floor < min_exp)
        {
            // as it turns out, taking log2 takes ~25 ns, but doing gsl_pow_int
            // takes 60 to 150 ns for powers above about 50
            // also, since m is in [0.5, 1), any cti > 200 or so will cause underflow)
            lp += cti * log2(q);
        }
        else
        {
            // this will have an exponent
            double m = frexp(q, &e);
            double em = gsl_pow_int(m, cti);
            s += (e * cti);
            p *= em;
            cti_floor -= cti;
        }
    }

    double lsm = log2(p) + s + lp;

    // assert(! isinf(ret));
    // assert(! isinf(ret2));
    return lsm / M_LOG2E;
}


// return log(P(x)), capping the value between -DBL_MAX and DBL_MAX
double ErrorEstimate::log_pdf_trunc(const double *x)
{
    double pri = gsl_ran_dirichlet_lnpdf(NUM_NUCS, this->prior_alpha, x);
    return this->log_likelihood(x) + (isinf(pri) ? isinf(pri) * DBL_MAX : pri);
}


//finds mode point of this posterior, returning true on success
//returns number of iterations
size_t ErrorEstimate::find_mode_point(double gradient_tolerance, 
                                      size_t max_iterations,
                                      const double *initial_point,
                                      bool *on_zero_boundary,
                                      bool verbose,
                                      double *mode_point) const
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

    gsl_multimin_fdfminimizer * first_minimizer =
        gsl_multimin_fdfminimizer_alloc(gsl_multimin_fdfminimizer_vector_bfgs,
                                        sphere_ndim);

    //shirk the initial point.  Just start in one of the four
    //quadrants that gets the highest value.
    double test_values[4];

    // precomputed sigmoid values representing 90%, 3.33%, 3.33%, 3.33% base composition.
    double test_r3_point[4][3] = {
        {2.6390567939010965, -0, 3.2958378660048289},
        {2.6390567939010965, -0, -3.2958378660048289},
        {-2.6390584010443274, 3.2958378660048289, -0},
        {-2.6390584010443274, -3.2958378660048289, -0}
    };
    
    for (size_t d = 0; d != 4; ++d)
    {
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



    //gsl_vector * numerical_gradient = gsl_vector_alloc(sphere_ndim);
    gsl_vector *last_point = gsl_vector_alloc(sphere_ndim);
    gsl_vector *point_delta = gsl_vector_alloc(sphere_ndim);

    gsl_vector_set_all(last_point, 0.0);

    gsl_vector *best_start_point = gsl_vector_alloc(sphere_ndim);
    gsl_vector_set(best_start_point, 0, test_r3_point[best_ind][0]);
    gsl_vector_set(best_start_point, 1, test_r3_point[best_ind][1]);
    gsl_vector_set(best_start_point, 2, test_r3_point[best_ind][2]);
    gsl_multimin_fdfminimizer_set(first_minimizer, &fdf_minimizer_function, best_start_point, 0.1, 0.1);
    gsl_vector_free(best_start_point);

    //double mode_gradient[4][3];

    //const double sigmoid_epsilon = 1e-10;

    Transformation::SigmoidVals sigmoid_vals[3];

    size_t mi;
    for (mi = 0; mi != max_iterations; ++mi)
    {

        //copy current point to last point
        gsl_vector_memcpy(last_point, first_minimizer->x);

        //take a step
        gsl_multimin_fdfminimizer_iterate(first_minimizer);

        gsl_vector *gradient = 
            gsl_multimin_fdfminimizer_gradient(first_minimizer);

        double x[3];
        x[0] = gsl_vector_get(first_minimizer->x, 0);
        x[1] = gsl_vector_get(first_minimizer->x, 1);
        x[2] = gsl_vector_get(first_minimizer->x, 2);
        
        Transformation::sigmoid_value_and_gradient(x, sigmoid_vals);

        //check if step is small enough
        gsl_vector_memcpy(point_delta, last_point);
        gsl_vector_sub(point_delta, first_minimizer->x);
        double last_step_size = gsl_blas_dnrm2(point_delta);

        if (verbose)
        {
            //print current stats
            
            Transformation::sigmoid_composition(sigmoid_vals, mode_point);
            printf(
                   "gslfdf1: iter: %04Zu, last_step_size: %10.8f, "
                   "log_neg: %20.15f  " 
                   "gradient: (%+8.6f, %+8.6f, %+8.6f)  "
                   "cur_mode: (%+8.6f, %+8.6f, %+8.6f, %+8.6f)  "
                   "r3_point: (%10.6f, %10.6f, %10.6f)\n"
                   ,
                   mi,
                   last_step_size,
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
        
        if (gsl_multimin_test_gradient(gradient, gradient_tolerance) == GSL_SUCCESS ||
            (mi > 20 && last_step_size < 1e-10)
            )
        {
            break;
        }

    }
    gsl_multimin_fdfminimizer_free(first_minimizer);
    gsl_vector_free(last_point);
    gsl_vector_free(point_delta);

    if (verbose)
    {
        printf("\n\n\n");
    }
    Transformation::sigmoid_composition(sigmoid_vals, mode_point);
    Transformation::boundary_point(sigmoid_vals, on_zero_boundary);

    return mi;
}
