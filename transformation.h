#include <gsl/gsl_vector.h>

class ErrorEstimate;

namespace Transformation
{

    int const VALUE = 1;
    int const GRADIENT = 2;
    
    struct PassingParams
    {
        const ErrorEstimate *error_estimate;
        double current_mode;
    };

    struct SigmoidVals
    {
        double p;
        double q;
        double pgrad;
        double qgrad;
        double pgrad_over_p;
        double qgrad_over_q;
        double log_p;
        double log_q;
    };


    void sigmoid_value_and_gradient(const double x[3],
                                    SigmoidVals sg[3]);

    void sigmoid_composition(SigmoidVals const sg[3],
                             double *comp);

    //determine whether posterior is concave along each of the four dimensions
    void boundary_point(const Transformation::SigmoidVals sg[3],
                        bool *on_zero_boundary);

    void sigmoid_gradient(const SigmoidVals sg[3],
                          double comp_gradient[4][3]);

    void sigmoid_log_dirichlet_gradient(const double alpha[4],
                                        const SigmoidVals sg[3],
                                        double *gradient);

    void composition_to_r3_sigmoid(const double *comp, double *r);

    double log_dirichlet(const double alpha[4],
                         const double x[4]);
    
    double log2_dirichlet(const double alpha[4], const double x[4]);

    double log_neg_posterior_value(const gsl_vector *r, void *params);

    void log_neg_posterior_gradient(const gsl_vector *r, 
                                    void *params, 
                                    gsl_vector *gradient);

    void log_neg_posterior_value_and_gradient(const gsl_vector *r, 
                                              void *params, 
                                              double *value, 
                                              gsl_vector *gradient);
    
/*     void log2_neg_gradient_numerical(const gsl_vector *r,  */
/*                                      void *params,  */
/*                                      double epsilon, */
/*                                      gsl_vector *numerical_gradient); */
    
};


