/*
// want a class that allows to sample, evaluate, memoize.
class AnalyticalIntegrand : public SamplingFunction
{
public:
    AnalyticalIntegrand(size_t ndim, bool _may_underflow) : 
        SamplingFunction(ndim, _may_underflow) { }

    virtual double marginal_cdf(double const xcoord, size_t const marg_dim) const = 0;
    virtual double inv_marginal_cdf(double const p, size_t const marg_dim) const = 0;
};


// implementation of multivariate gaussian
class Gaussian : public AnalyticalIntegrand
{

    gsl_vector * mean, * cdf_lowbound, * cdf_highbound;
    gsl_matrix * covariance, * covariance_inv, * cholesky, * L;

    // size_t ndim;
    
    double factor;
    double normalization_constant;

    gsl_rng * gsl_rand_gen;
    gsl_error_handler_t * default_handler;

    bool initialized;

public:
    
    Gaussian(double const* _mean, double const* _covariance, size_t const _ndim,
             bool _may_underflow);

    ~Gaussian();
    void Init();

    void get_mean(double * _mean);
    void set_mean(double const* _mean);

    void get_covariance(double * _covariance);
    void set_covariance(double const* _covariance);

    void update_mean(double const* _mean);
    void update_covariance(double const* _covariance);

    double pdf(double const* x);
    double log_pdf(double const* x);

    double marginal_cdf(double const xcoord, size_t const marg_dim) const;
    // double marginal_truncated_cdf(double const xcoord, size_t const marg_dim);
    double inv_marginal_cdf(double const p, size_t const marg_dim) const;
    // double inv_marginal_truncated_cdf(double const p, size_t const marg_dim);
    void sample(double * x) const;
    void sample_conditioned(double const* x, double * y) { }
    void print_params(FILE * fh) const;
};
*/
