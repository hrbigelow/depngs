#include "transformation.h"
#include "stats_tools.h"
#include "error_estimate.h"

#include <cassert>

#include <gsl/gsl_sf_trig.h>
#include <gsl/gsl_sf_exp.h>
#include <gsl/gsl_sf_log.h>
#include <gsl/gsl_complex_math.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_randist.h>

namespace Transformation
{


    const double sigmoid_truncation_epsilon = 1e-10;
    const double sigmoid_truncation_epsilon_inv = 1.0 - sigmoid_truncation_epsilon;
    const double log_epsilon = gsl_sf_log(sigmoid_truncation_epsilon);
    const double log_epsilon_inv = gsl_sf_log(sigmoid_truncation_epsilon_inv);

    /*
      filipe's iterative transformation: take a pair of additive inverses in one variable, and
      iteratively add another variable by multiplying by the sum of a pair:
      P(x) + Q(x) = 1
      P(x) + (P(y)+Q(y))Q(x) = 1
      (P(z) + Q(z))P(x) + (P(y)+Q(y))Q(x) = 1

      P(z)P(x) + Q(z)P(x) + P(y)Q(x) + Q(y)Q(x) = 1

      Let 
      P(x) + 1/(1 + exp(-x)), 
      Q(x) = exp(-x) / (1 + exp(-x))

      
    */
    void sigmoid_value_and_gradient(const double x[3],
                                    Transformation::SigmoidVals sg[3])
    {
        double enegx;

        const double epsilon = sigmoid_truncation_epsilon;
        const double epsilon_inv = 1.0 - epsilon;

        //here we assert that if the value p[d] is near 0 or 1,
        //the partial derivatives p'[d] and q'[d] we will deem
        //them to be zero, as if we'd already reached infinity
        //analytically, the sg and dirichlet cancel, and the
        //partial approaches a positive constant towards infinity.
        //but, for purposes of handling the dirichlet prior, we
        //make this correction.
        Transformation::SigmoidVals *sgp;
        
        for (size_t d = 0; d != 3; ++d)
        {
            // gsl_error_handler_t * original_handler = gsl_set_error_handler_off();
            enegx = gsl_sf_exp(- x[d]);
            // gsl_set_error_handler(original_handler);

            sgp = & sg[d];

            sgp->p = 1.0 / (1.0 + enegx);

            if (sgp->p < epsilon)
            {
                sgp->p = epsilon;
                sgp->pgrad = 0.0;
                sgp->pgrad_over_p = 0.0;
                sgp->log_p = log_epsilon;

                sgp->q = epsilon_inv;
                sgp->qgrad = 0.0;
                sgp->qgrad_over_q = 0.0;
                sgp->log_q = log_epsilon_inv;
            }
            else if (sgp->p > epsilon_inv)
            {
                sgp->p = epsilon_inv;
                sgp->pgrad = 0.0;
                sgp->pgrad_over_p = 0.0;
                sgp->log_p = log_epsilon_inv;

                sgp->q = epsilon;
                sgp->qgrad = 0.0;
                sgp->qgrad_over_q = 0.0;
                sgp->log_q = log_epsilon;
            }
            else 
            {
                sgp->q = 1.0 - sgp->p;
                sgp->pgrad = (gsl_isinf(enegx) == 1) ? 0.0 : enegx / gsl_pow_2(1.0 + enegx);
                // sgp->pgrad = gsl_ran_logistic_pdf(x[d], 1.0);
                sgp->qgrad = - sgp->pgrad;
                
                double log_1plus_enegx = gsl_sf_log_1plusx(enegx);
                sgp->log_p = -log_1plus_enegx;
                sgp->log_q = - x[d] - log_1plus_enegx;
                
                sgp->pgrad_over_p = sgp->q;
                sgp->qgrad_over_q = - sgp->p;
            }

            assert(! isnan(sgp->p));
            assert(! isnan(sgp->q));
            assert(! isnan(sgp->pgrad_over_p));
            assert(! isnan(sgp->qgrad_over_q));
            assert(! isnan(sgp->pgrad));
            assert(! isnan(sgp->qgrad));
            assert(! isnan(sgp->log_p));
            assert(! isnan(sgp->log_q));
            
        }

    }                     



    // compute the original 4D point from the sigmoid transformed point
    void sigmoid_composition(const Transformation::SigmoidVals sg[3],
                             double *comp)
    {
        comp[0] = sg[2].p * sg[0].p;
        comp[1] = sg[2].q * sg[0].p;
        comp[2] = sg[1].p * sg[0].q;
        comp[3] = sg[1].q * sg[0].q;

        // !!! might be a mistake to comment this out...
        // if (! (normalized(comp, 4, 1e-10) && all_positive(comp, 4)))
        // {
        //     fprintf(stderr, "r3_to_composition_sigmoid: invalid composition\n");
        //     exit(2);
        // }
        
    }

    //if the transformation has reached a point where the gradients are flat,
    //this is deemed to be 'infinity', which means the point represents the boundary
    void boundary_point(const Transformation::SigmoidVals sg[3],
                        bool *on_zero_boundary)
    {
        const double& e = Transformation::sigmoid_truncation_epsilon;

        on_zero_boundary[0] = sg[2].p == e || sg[0].p == e;
        on_zero_boundary[1] = sg[2].q == e || sg[0].p == e;
        on_zero_boundary[2] = sg[1].p == e || sg[0].q == e;
        on_zero_boundary[3] = sg[1].q == e || sg[0].q == e;

    }

    void sigmoid_gradient(const Transformation::SigmoidVals sg[3],
                          double comp_gradient[4][3])
    {
        comp_gradient[0][0] = sg[2].p * sg[0].pgrad;
        comp_gradient[0][1] = 0.0;
        comp_gradient[0][2] = sg[0].p * sg[2].pgrad;
        
        comp_gradient[1][0] = sg[2].q * sg[0].pgrad;
        comp_gradient[1][1] = 0.0;
        comp_gradient[1][2] = sg[0].p * sg[2].qgrad;

        comp_gradient[2][0] = sg[1].p * sg[0].qgrad;
        comp_gradient[2][1] = sg[0].q * sg[1].pgrad;
        comp_gradient[2][2] = 0.0;

        comp_gradient[3][0] = sg[1].q * sg[0].qgrad;
        comp_gradient[3][1] = sg[0].q * sg[1].qgrad;
        comp_gradient[3][2] = 0.0;

    }


    //compute log(~Dir(sigmoid(x)))
    double sigmoid_log_dirichlet(const double alpha[4],
                                 const Transformation::SigmoidVals sg[3])
    {
        double i[4];
        for (size_t c = 0; c != 4; ++c)
        {
            i[c] = alpha[c] - 1.0;
        }
        double value =
            + i[0] * (sg[2].log_p + sg[0].log_p)
            + i[1] * (sg[2].log_q + sg[0].log_p)
            + i[2] * (sg[1].log_p + sg[0].log_q)
            + i[3] * (sg[1].log_q + sg[0].log_q);

        return value;
    }
    
    //compute d/dx of log(~Dir(sigmoid(x)))
    void sigmoid_log_dirichlet_gradient(const double alpha[4],
                                        const Transformation::SigmoidVals sg[3],
                                        double *gradient)
    {
        double i[4];
        for (size_t c = 0; c != 4; ++c)
        {
            i[c] = alpha[c] - 1.0;
        }
        gradient[0] = 
            + i[0] * sg[0].pgrad_over_p 
            + i[1] * sg[0].pgrad_over_p 
            + i[2] * sg[0].qgrad_over_q 
            + i[3] * sg[0].qgrad_over_q;
        
        gradient[1] = i[2] * sg[1].pgrad_over_p + i[3] * sg[1].qgrad_over_q; 
        gradient[2] = i[0] * sg[2].pgrad_over_p + i[1] * sg[2].qgrad_over_q;
    }

    /*
      Let 

      Pinv(x) = -log(1/y - 1)

      c0 = P(z)P(x)
      c1 = Q(z)P(x)
      c2 = P(y)Q(x)
      c3 = Q(y)Q(x)

      so:

      x = Pinv(c0 + c1)
      y = Pinv(c2 / Q(x)) = Pinv(c2 / (c2 + c3)) = -log(1 + c3/c2 - 1)
      z = Pinv(c0 / P(x)) = Pinv(c0 / (c0 + c1)) = -log(1 + c1/c0 - 1)
      
    */

    //transform normalized 4D coordinates to r3 using a sigmoid
    //this needs to apply the truncation logic 
    //gsl_sf_log returns -nan if zero or close to zero with status 1,
    //inf if argument is larger than representable, with status 0
    /*
    void composition_to_r3_sigmoid(const double* c, double *r)
    {
        
        // gsl_error_handler_t *original_handler = gsl_set_error_handler_off();
        double arg[3];
        arg[0] = 1.0 / (c[0] + c[1]) - 1.0;
        arg[1] = (c[3] == 0 && c[2] == 0) ? 1 : c[3] / c[2];
        arg[2] = (c[1] == 0 && c[0] == 0) ? 1 : c[1] / c[0];

        int status;
        gsl_sf_result result;
        status = gsl_sf_log_e(arg[0], & result);
        r[0] = status ? FLT_MAX : ((gsl_isinf(result.val) == 1) ? -FLT_MAX : -result.val);

        status = gsl_sf_log_e(arg[1], & result);
        r[1] = status ? FLT_MAX : ((gsl_isinf(result.val) == 1) ? -FLT_MAX : -result.val);

        status = gsl_sf_log_e(arg[2], & result);
        r[2] = status ? FLT_MAX : ((gsl_isinf(result.val) == 1) ? -FLT_MAX : -result.val);

        // gsl_set_error_handler(original_handler);
    }
    */

    // may be NaN, this is okay
    /*
    double log_dirichlet(const double *alpha, const double *x)
    {
        // double logs[] = { 
        //     gsl_sf_log(x[0]),
        //     gsl_sf_log(x[1]),
        //     gsl_sf_log(x[2]),
        //     gsl_sf_log(x[3])
        // };
        
        
        double ld =
            + (alpha[0] == 1.0 ? 0 : (alpha[0] - 1.0) * gsl_sf_log(x[0]))
            + (alpha[1] == 1.0 ? 0 : (alpha[1] - 1.0) * gsl_sf_log(x[1]))
            + (alpha[2] == 1.0 ? 0 : (alpha[2] - 1.0) * gsl_sf_log(x[2]))
            + (alpha[3] == 1.0 ? 0 : (alpha[3] - 1.0) * gsl_sf_log(x[3]));

        return ld;
    }
    */


    /*
    double log2_dirichlet(const double alpha[4], const double x[4])
    {
        double ld =
            + (alpha[0] - 1.0) * log2(x[0])
            + (alpha[1] - 1.0) * log2(x[1])
            + (alpha[2] - 1.0) * log2(x[2])
            + (alpha[3] - 1.0) * log2(x[3]);
        
        // if (ld != 0)
        // {
        //     assert(ld / ld2 < 1.00001);
        //     assert(ld2 / ld < 1.00001);
        // }
        return ld;
        
    }
    */


    

    void log_neg_posterior_aux(const gsl_vector *r, 
                               int eval_type, 
                               double *neg_log_value,
                               gsl_vector *neg_gradient_vec,
                               void *params)
    {
        double comp[4];
        double comp_gradient[4][3];
        double prior_gradient[3];
        double likelihood_gradient[4];
        double gradient[3];

        double log_prior;
        double log_likelihood;

        double x[3];
        x[0] = gsl_vector_get(r, 0);
        x[1] = gsl_vector_get(r, 1);
        x[2] = gsl_vector_get(r, 2);

        Transformation::PassingParams *pp = 
            static_cast<Transformation::PassingParams *>(params);

        Transformation::SigmoidVals sigmoid_vals[3];

        //smallest allowed distance from 0 or 1.  any closer to this,
        //and the sigmoid is considered to be equal to 0 or 1 with
        //derivative equal to 0

        sigmoid_value_and_gradient(x, sigmoid_vals);

        sigmoid_composition(sigmoid_vals, comp);

        if (eval_type & Transformation::VALUE)
         //if (1)
        {
            log_likelihood = pp->error_estimate->log_likelihood(comp) - pp->current_mode;
            log_prior = 
                pp->error_estimate->uniform_prior
                ? 0
                : sigmoid_log_dirichlet(pp->error_estimate->prior_alpha,
                                        sigmoid_vals);

            *neg_log_value = -1.0 * (log_likelihood + log_prior);
        }

        if (eval_type & Transformation::GRADIENT)
        // if (1)
        {
            sigmoid_gradient(sigmoid_vals, comp_gradient);

            if (pp->error_estimate->uniform_prior)
            {
                std::fill(prior_gradient, prior_gradient + 3, 0.0);
            }
            else
            {
                sigmoid_log_dirichlet_gradient(pp->error_estimate->prior_alpha,
                                               sigmoid_vals, prior_gradient);
            }
            
            pp->error_estimate->log_likelihood_gradient(comp, likelihood_gradient);

            for (size_t ri = 0; ri != 3; ++ri)
            {
                gradient[ri] = 0.0;
                for (size_t c = 0; c != 4; ++c)
                {
                    gradient[ri] += likelihood_gradient[c] * comp_gradient[c][ri];
                }
                gradient[ri] += prior_gradient[ri];
                gsl_vector_set(neg_gradient_vec, ri, gradient[ri] * -1.0);
            }
        }
        //if (1)
        if (0)
        {
            printf("etype: %i, "
                   "val: %10.8f, ll: %10.8f, prior: %10.8f, "
                   "r3: %10.8f, %10.8f, %10.8f, "
                   "prg: %10.8f, %10.8f, %10.8f, "
                   "grd: %10.8f, %10.8f, %10.8f, "
                   "composition: %10.8f, %10.8f, %10.8f, %10.8f\n",
                   eval_type,
                   *neg_log_value,
                   log_likelihood, log_prior,
                   x[0], x[1], x[2],
                   prior_gradient[0], prior_gradient[1], prior_gradient[2],
                   gradient[0], gradient[1], gradient[2],
                   comp[0], comp[1], comp[2], comp[3]);
        }
    }

    double log_neg_posterior_value(const gsl_vector *r, 
                                   void *params)
    {
        double neg_log_value;
        gsl_vector *neg_log_gradient = gsl_vector_alloc(3);

        log_neg_posterior_aux(r, Transformation::VALUE, &neg_log_value, 
                              neg_log_gradient, params);

        gsl_vector_free(neg_log_gradient);

        return neg_log_value;
    }


    void log_neg_posterior_gradient(const gsl_vector *r, 
                                    void *params, 
                                    gsl_vector *neg_gradient)
    {
        double neg_log_value;
        log_neg_posterior_aux(r, Transformation::GRADIENT, 
                              &neg_log_value, neg_gradient, params);
    }

    void log_neg_posterior_value_and_gradient(const gsl_vector *r, 
                                              void *params, 
                                              double *neg_log_value, 
                                              gsl_vector *neg_gradient)
    {
        log_neg_posterior_aux(r, Transformation::VALUE | Transformation::GRADIENT, 
                              neg_log_value, neg_gradient, params);
    }

} // END namespace Transformation
