//Produce the next sample point from the posterior
//!!! need a way to get the mode point, pack it, and calculate
//the log2_scaling_factor for that point.
std::vector<std::pair<BaseComp, REAL> > 
SliceSampling(ErrorEstimate const& ee, 
              BASE_COMP_COUNTS const& base_comp_counts,
              BaseComp const& starting_point,
              REAL log2_scaling_factor,
              int const initial_range,
              size_t const num_samples)
{
    size_t const ndim = 3;
    size_t const bits = sizeof(double) * 8;

    //xprime and yprime are the proposed new point for slice sampling
    //they are only accepted if yprime < y.  once accepted, the slice
    //sampling continues, using (xprime, yprime) as the new (x, y)

    REAL yprime, y;
    mpz_t U, N, B, x, xprime;
    mpf_t uniform;

    mpz_init(U);
    mpz_init(N);
    mpz_init2(B, 0);
    mpz_init(xprime);

    mpf_init(uniform);
    
    mpz_setbit(B, ndim * bits);

//     mpf_out_str(stderr, 10, 50, Bf);
//     fprintf(stderr, "\n");

    BaseComp sample_composition;

    gmp_randstate_t rand_state;
    gmp_randinit_default(rand_state);

    mpf_urandomb(uniform, rand_state, sizeof(REAL) * 8);

    y = ee.ScaledPosterior(base_comp_counts, starting_point, 
                           log2_scaling_factor) * mpf_get_d(uniform);

    PackBaseComposition(ee, starting_point, & x);

    sample_composition = UnpackBaseComposition(ee, x);

    int current_range;

    std::vector<std::pair<BaseComp, REAL> > sample_points;
    sample_points.resize(num_samples);

    //generate the sample points
    for (size_t si = 0; si != num_samples; ++si)
    {
        
        mpz_urandomb(U, rand_state, ndim * bits);
        current_range = initial_range;
        do
        {
            //goes through the space of integers
            mpz_urandomb(N, rand_state, current_range);

//             fprintf(stderr, "= x: ");
//             mpz_out_str(stderr, 2, x);
//             fprintf(stderr, "\n");
            
            mpz_sub(xprime, x, U);

//             fprintf(stderr, "- U: ");
//             mpz_out_str(stderr, 2, xprime);
//             fprintf(stderr, "\n");

            mpz_mod(xprime, xprime, B);

//             fprintf(stderr, "M B: ");
//             mpz_out_str(stderr, 2, xprime);
//             fprintf(stderr, "\n");

            mpz_xor(xprime, xprime, N);

//             fprintf(stderr, "^ N: ");
//             mpz_out_str(stderr, 2, xprime);
//             fprintf(stderr, "\n");

            mpz_add(xprime, xprime, U);

//             fprintf(stderr, "+ U: ");
//             mpz_out_str(stderr, 2, xprime);
//             fprintf(stderr, "\n");

            mpz_mod(xprime, xprime, B);

//             fprintf(stderr, "M B: ");
//             mpz_out_str(stderr, 2, xprime);
//             fprintf(stderr, "\n");

//             fflush(stderr);

            --current_range;

            sample_composition = UnpackBaseComposition(ee, xprime);

//             mpz_t testpack;
//             PackBaseComposition(ee, sample_composition, & testpack);
//             BaseComp testcomp = UnpackBaseComposition(ee, testpack);


            yprime = ee.ScaledPosterior(base_comp_counts, sample_composition, 
                                        log2_scaling_factor);
            
            int xxx = 0;
        }
        while ((mpz_cmp(xprime, x) != 0) && yprime < y);

        sample_points[si] = std::make_pair(sample_composition, yprime);

        //adopt current
        mpz_set(x, xprime);

        mpf_urandomb(uniform, rand_state, sizeof(REAL) * 8);
        y = yprime * mpf_get_d(uniform);
        
    }
    return sample_points;
}








//Convert a basecomp into a packed coordinate
void PackBaseComposition(ErrorEstimate const& ee, 
                         BaseComp const& basecomp,
                         mpz_t * packed)
{

    size_t const ndim = 3;
    size_t const bits = sizeof(double) * 8;

    mpz_t coord[ndim];
    mpf_t coordf[ndim];

    //expand the coordinates to the hypercube
    double cube[ndim], slice[ndim];

    slice[0] = basecomp.A;
    slice[1] = basecomp.C;
    slice[2] = basecomp.G;

    ee.expand_to_hypercube(slice, cube);

    //rescale the number from [0,1] to [0, INT_MAX]
    uint64_t intcube[ndim];
    for (size_t i = 0; i != ndim; ++i)
    {
        intcube[i] = static_cast<uint64_t>(cube[i] * ~ static_cast<uint64_t>(0));
    }
    Hilbert_to_int(intcube, ndim, 64, *packed);

}


BaseComp UnpackBaseComposition(ErrorEstimate const& ee, 
                               mpz_t const& packed)
{
    
    const size_t ndim = 3;

    uint64_t coord[ndim];

    //expand the coordinates to the hypercube
    double cube[ndim], slice[ndim];

    int_to_Hilbert(packed, 64, coord, ndim);

    for (size_t i = 0; i != ndim; ++i)
    {
        cube[i] = static_cast<double>(coord[i]) / static_cast<double>(~static_cast<uint64_t>(0));
    }
    
    ee.contract_from_hypercube(cube, slice);
    
    return BaseComp(slice[0], slice[1], slice[2],
                    1.0 - slice[0] - slice[1] - slice[2]);
    
}






        //Find the mode point
        gsl_multimin_function minimizer_function;

        //Bounds full_bounds[] = { { 0.0, 1.0 }, { 0.0, 1.0 }, { 0.0, 1.0 } };

        minimizer_function.n = 3;
        minimizer_function.f = &gsl_simplex_hypercube_posterior;
        //minimizer_function.f = &gsl_simplex_posterior;

        
        bool add_teepee_component = true; //adds a very gently sloping teepee to the posterior to help find the peak

        const int initial_call_count = 0;

        PosteriorParams posterior_params = {
            &error_estimate, &base_comp_counts, 0.0, 
            static_cast<std::vector<WeightedSample> *>(NULL), 
            { { 0.0, 1.0 }, { 0.0, 1.0 }, { 0.0, 1.0 } },
            scale_increment,
            scale_increment,
            initial_call_count,
            0.0, //initial mode point
            { 0.0, 0.0, 0.0 } // mode point placeholder
        };
        

        minimizer_function.params = &posterior_params;
        gsl_vector *ss, *x, *x_cur, *x_prev, *x_jump;
        gsl_multimin_fminimizer * minimizer = 
            gsl_multimin_fminimizer_alloc(gsl_multimin_fminimizer_nmsimplex2, 3);

        x = gsl_vector_alloc(3);
        gsl_vector_set(x, 0, 0.25);
        gsl_vector_set(x, 1, 0.25);
        gsl_vector_set(x, 2, 0.25);

        ss = gsl_vector_alloc(3);
        gsl_vector_set_all(ss, 0.1);

        x_cur = gsl_vector_alloc(3);
        x_prev = gsl_vector_alloc(3);
        x_jump = gsl_vector_alloc(3);

        gsl_multimin_fminimizer_set(minimizer, &minimizer_function, x, ss);

       
        double prev_simplex_size = 10.0;
        double simplex_size = 1.0;
        double jump_size = 1.0;

        int max_multimin_iterations = 100000;
        int multimin_iteration = 0;
        while (gsl_multimin_test_size(simplex_size, mode_tolerance) == GSL_CONTINUE
               && multimin_iteration < max_multimin_iterations)
        {
            prev_simplex_size = simplex_size;
            gsl_multimin_fminimizer_iterate(minimizer);
            simplex_size = gsl_multimin_fminimizer_size(minimizer);
            ++multimin_iteration;
            if (multimin_iteration % 10000 == 0)
            {
                //                 printf("simplex_size: %20.20g, iteration: %i\n",
                //                        simplex_size, static_cast<int>(multimin_iteration));
            }
        }

        //         printf("Found mode in %i iterations with simplex size %20.20f\n",
        //                multimin_iteration, simplex_size);

        double mode[] = {
            gsl_vector_get(minimizer->x,0),
            gsl_vector_get(minimizer->x,1),
            gsl_vector_get(minimizer->x,2)
        };


        double expanded_mode[] = {
            gsl_vector_get(minimizer->x,0),
            gsl_vector_get(minimizer->x,1),
            gsl_vector_get(minimizer->x,2)
        };

        //error_estimate.expand_to_hypercube(mode, expanded_mode);
        error_estimate.contract_from_hypercube(expanded_mode, mode);

        //error_estimate.contract_from_hypercube(expanded_mode, mode);

        BaseComp mode_comp(mode[0], mode[1], mode[2], 1.0 - mode[0] - mode[1] - mode[2]);
        
        //deallocate minimizer and x vector
        gsl_multimin_fminimizer_free(minimizer);
        gsl_vector_free(x);
        gsl_vector_free(ss);

        posterior_params.best_log2_mode = error_estimate.Log2Posterior(base_comp_counts, mode_comp);

        //initially set the shift factor to the best log2 mode.
        posterior_params.shift_factor = posterior_params.best_log2_mode;

        //         printf("Mode point: %f %f %f %f, Log2Posterior at mode = %f\n",
        //                mode_comp.A, mode_comp.C, mode_comp.G, mode_comp.T, posterior_params.shift_factor);

        posterior_params.mode_point[0] = mode[0];
        posterior_params.mode_point[1] = mode[1];
        posterior_params.mode_point[2] = mode[2];

        

        
//accepts a point in the unit hypercube, condenses it,
//evaluates the sensored log2 posterior (0 for any underflow values)
//in the legal region
void cuba_posterior(int const* ndim, double const cube_coord[],
                    int const* ncomp, void * posterior_params, 
                    double output[], double * weight)
{
    
    double contraction_scale = 1.0 / 6.0;
    PosteriorParams * pp = static_cast<PosteriorParams *>(posterior_params);
    double pyramid_coord[3];
    pp->ee->contract_from_hypercube(cube_coord, pyramid_coord);

    BaseComp composition(pyramid_coord[0], pyramid_coord[1], pyramid_coord[2], 
                         1.0 - pyramid_coord[0] - pyramid_coord[1] - pyramid_coord[2]);

    REAL log2_posterior = pp->ee->Log2Posterior(composition);
        
    REAL log2_shifted = log2_posterior - pp->shift_factor;
    
    REAL log2_scaled = log2_shifted * pp->scale_factor;

    REAL ceiling = std::numeric_limits<REAL>::max_exponent;
    if (log2_scaled > ceiling)
    {
        log2_scaled = ceiling - 10.0;
        //taper + log2f((log2_scaled - taper) + 1);
    }

    //log2_scaled needs to be exponentiated.  if it is out of range either high or low,
    //it will be truncated.
    REAL safe_log2posterior = exp2l(log2_scaled);
   
    
    // !!! warning: not using contraction scale
    output[0] = safe_log2posterior;

    assert(output[0] != std::numeric_limits<REAL>::infinity() &&
           output[0] != std::numeric_limits<REAL>::quiet_NaN() &&
           output[0] != std::numeric_limits<REAL>::signaling_NaN());


    if (log2_posterior > pp->best_log2_mode)
    {
        pp->best_log2_mode = log2_posterior;
        pp->mode_point[0] = pyramid_coord[0];
        pp->mode_point[1] = pyramid_coord[1];
        pp->mode_point[2] = pyramid_coord[2];

        //printf("%f\n", log2_shifted);
        //fflush(stdout);
        //printf("cuba_posterior above bound: %f\n", output[0]);
    }

    //assert(output[0] <= 1.0);
    if (pp->call_count % static_cast<int>(1e6 * pp->scale_increment) == 0)
    {
        pp->scale_factor += pp->scale_increment;
        pp->scale_factor = std::min(pp->scale_factor, 1.0);
        //printf("scale_factor: %f\n", pp->scale_factor);
    }
    ++pp->call_count;

    if (pp->samples != NULL)
    {
        (*pp->samples).push_back(WeightedSample(4, pyramid_coord, output[0], *weight));
    }
}





//simple wrapper for the cuba posterior to work with the divonne function
//phase is not used, and weight
void divonne_posterior(int const* ndim, double const x[],
                       int const* ncomp, void * posterior_params, 
                       double output[], int const* phase)
{
    double weight = 1.0;
    return cuba_posterior(ndim, x, ncomp, posterior_params, output, &weight);
}



void cuhre_posterior(int const* ndim, double const x[],
                     int const* ncomp, void * posterior_params, 
                     double output[], double * dummy)
{
    double weight = 1.0;
    return cuba_posterior(ndim, x, ncomp, posterior_params, output, &weight);
}



void cuba_posterior(int const* ndim, double const x[],
                    int const* ncomp, void * posterior_params, 
                    double output[], double * weight);


void divonne_posterior(int const* ndim, double const x[],
                       int const* ncomp, void * posterior_params, 
                       double output[], int const* phase);


void cuhre_posterior(int const* ndim, double const x[],
                     int const* ncomp, void * posterior_params, 
                     double output[], double * dummy);
                        



void PrintCDFComparison(FILE * out_fh, 
                        Gaussian const* gaussian,
                        WEIGHTED_SAMPLE_MAP & samples,
                        double const* quantiles, 
                        size_t const num_quantiles,
                        size_t const num_dimensions)
{

    std::vector<WeightedSample *> samples_vec;
    WEIGHTED_SAMPLE_MAP::iterator sample_iter;

    for (sample_iter = samples.begin();
         sample_iter != samples.end(); ++sample_iter)
    {
        WeightedSample * ws = & (*sample_iter).second;
        samples_vec.push_back(ws);
    }
    

    //print out a number of quantiles
    for (size_t d = 0; d != num_dimensions; ++d)
    {
        MarginalCumulativeDistribution(d, &samples_vec);
    }

    for (size_t q = 0; q != num_quantiles; ++q)
    {
        bool is_lower_bound = quantiles[q] < 0.5;

        for (size_t d = 0; d != num_dimensions; ++d)
        {
            
            REAL inverse_cdf_est = 
                FindIntegralBound(&samples_vec, d,
                                  quantiles[q], 
                                  is_lower_bound);
            
            double inverse_cdf = gaussian->inv_marginal_cdf(quantiles[q], d);
            double quantile_est = gaussian->marginal_cdf(inverse_cdf_est, d);

            fprintf(stdout, "\t%10.8f\t%10.8f\t%10.8f\t%10.8f", 
                    quantiles[q],
                    quantile_est - quantiles[q],
                    inverse_cdf,
                    inverse_cdf_est - inverse_cdf);
        }
        fprintf(stdout, "\n");

    }
}


void find_integral_bounds(std::vector<double *> * points,
                          size_t sort_dimension,
                          double const* quantiles,
                          size_t num_quantiles,
                          double * quantile_values)
{
    size_t num_points = (*points).size();
    std::sort((*points).begin(), (*points).end(), SortDimension(sort_dimension));
    for (size_t f = 0; f != num_quantiles; ++f)
    {
        size_t cut_point;
        if (quantiles[f] < 0.5)
        {
            cut_point = 
                std::floor(quantiles[f] * static_cast<double>(num_points));
        }
        else
        {
            cut_point = 
                std::ceil(quantiles[f] * static_cast<double>(num_points));
        }
            
        quantile_values[f] = (*points)[cut_point - 1][sort_dimension];
    }
}



ProposalGaussian::ProposalGaussian(Integrand * _integrand,
                                   double const* _mu,
                                   double const* _sigma, 
                                   size_t const _ndim) : 
    integrand(_integrand), ndim(_ndim)
{ 
    this->mu = new double[ndim];
    this->sigma = new double[ndim*ndim];
//     this->z_tau = new double[ndim];
//     this->z_star = new double[ndim];
    std::copy(_mu, _mu + ndim, this->mu);
    std::copy(_sigma, _sigma + (ndim * ndim), this->sigma);
    gsl_rand_gen = gsl_rng_alloc(gsl_rng_taus);
    //gsl_rng_set(const gsl_rng * r, unsigned long int s)
}

ProposalGaussian::~ProposalGaussian()
{
    if (this->sigma != NULL)
    {
        delete this->sigma;
        this->sigma = NULL;
    }
    if (this->mu != NULL)
    {
        delete this->mu;
        this->mu = NULL;
    }
//     if (this->z_star != NULL)
//     {
//         delete this->z_star;
//         this->z_star = NULL;
//     }
    gsl_rng_free(gsl_rand_gen);
}


void ProposalGaussian::propose(double const* z_tau, 
                               double * z_star)
{
    for (size_t d = 0; d != this->ndim; ++d)
    {
        z_star[d] = 
            gsl_ran_gaussian(this->gsl_rand_gen, this->sigma[d]) + z_tau[d];
    }
}


double ProposalGaussian::accept(double const* z_star,
                                double const* z_tau)
{
    double y_star = (*this->integrand)(z_star, this->ndim);
    double y_tau = (*this->integrand)(z_tau, this->ndim);
    
    if (y_tau > 0.0)
    {
        return std::min(1.0, y_star / y_tau);
    }
    else
    {
        return 0.0;
    }
    
}



class ProposalGaussian : public ProposalDistribution
{

    Integrand * integrand;
    size_t ndim;
    double * mu;
    double * sigma;

    gsl_rng * gsl_rand_gen;

 public:
    ProposalGaussian(Integrand * _integrand,
                     double const* _mu,
                     double const* _sigma, size_t const _ndim);

    ~ProposalGaussian();

    void propose(double const* z_tau, 
                 double * z_star);

    double accept(double const* z_star, 
                  double const* z_tau);
};



/*
  
 */


class MetropolisStepping
{


 public:
    virtual void propose(double const* z_tau, 
                         double * z_star) = 0;

    virtual double accept(double const* z_star, 
                          double const* z_tau) = 0;

    virtual void get_mean(double * mean) = 0;
    virtual void update_mean(double const* mean) = 0;

    virtual void get_covariance(double * covariance) = 0;
    virtual void update_covariance(double const* covariance) = 0;

};



//wraps the proposal distribution and integrand in a framework for
//proposing, accepting, and updating.
class IndependenceChainMH : public MetropolisStepping
{

    size_t ndim;
    bool log2_integrand;

 public:


    IndependenceChainMH(Integrand * _integrand,
                       double const* _mu,
                       double const* _sigma, 
                       size_t const _ndim,
                       bool _log2_integrand);

    ~IndependenceChainMH();

    bool in_bounds(double const* z_star);


    void get_mean(double * mean);
    void update_mean(double const* mean);

    void get_covariance(double * covariance);
    void update_covariance(double const* covariance);

    void Init();
    
};



IndependenceChainMH::IndependenceChainMH(Integrand * _integrand,
                                         SamplingFunction * _proposal,
                                         double const* _mu,
                                         double const* _sigma, 
                                         size_t const _ndim,
                                         bool _log2_integrand) : 
    integrand(_integrand), proposal(_proposal), 
    ndim(_ndim), log2_integrand(_log2_integrand)
{ 
}

IndependenceChainMH::~IndependenceChainMH()
{
}



bool IndependenceChainMH::in_bounds(double const* z_star)
{
    double sum = 0.0;
    bool in_bounds = true;
    for (size_t d = 0; d != this->ndim; ++d)
    {
        sum += z_star[d];
        in_bounds = in_bounds && z_star[d] >= 0.0;
    }
    in_bounds = in_bounds && sum <= 1.0;
    return in_bounds;
}



void IndependenceChainMH::get_mean(double * mean)
{
    this->proposal.get_mean(mean);
}


void IndependenceChainMH::update_mean(double const* mean)
{
    this->proposal.update_mean(mean);
}


void IndependenceChainMH::get_covariance(double * covariance)
{
    this->proposal.get_covariance(covariance);
}


void IndependenceChainMH::update_covariance(double const* covariance)
{
    this->proposal.update_covariance(covariance);
}



void IndependenceChainMH::Init()
{
    this->proposal.Init();
}




/*
void ErrorEstimate::find_mode_point2(double const* startx, double * found_mode)
{
    
    const int ndim = 4;
    int n_output_dim = this->base_comp_counts.size();

    std::copy(startx, startx + ndim, found_mode);

    double lb[] = { 0.0, 0.0, 0.0, 0.0 };
    double ub[] = { 1.0, 1.0, 1.0, 1.0 };
    //double ub[] = { 0.35, 0.35, 0.35, 0.35 };
    double A[] = { 1.0, 1.0, 1.0, 1.0 };
    int k = 1;
    double b[] = { 1.0 };
//     double * weights;
    int max_iterations = 100000;
//     double *options;
//     double *info;
    double * target_value = new double[n_output_dim];
    std::fill(target_value, target_value + n_output_dim, 0.0);

    dlevmar_blec_dif(&levmar_log2posterior, found_mode, target_value, ndim, n_output_dim,
                     lb, ub,
                     A, b, k,
                     NULL,
                     max_iterations, 
                     NULL,
                     NULL, NULL, NULL, this);

    delete target_value;

}
*/



/*
    fprintf(stdout, "mode after crude method: %g\t%g\t%g\n",
            gsl_vector_get(minimizer->x, 0),
            gsl_vector_get(minimizer->x, 1),
            gsl_vector_get(minimizer->x, 2));

    double step_size = 1e-2;

    gsl_multimin_fdfminimizer * fdf_minimizer =
        gsl_multimin_fdfminimizer_alloc(gsl_multimin_fdfminimizer_vector_bfgs2, ndim);
    
    gsl_multimin_function_fdf fdf_minimizer_function;

    fdf_minimizer_function.n = ndim;
    fdf_minimizer_function.params = & params;
    fdf_minimizer_function.f = &gsl_simplex_log2posterior;
    fdf_minimizer_function.df = &gsl_simplex_log2posterior_gradient;
    fdf_minimizer_function.fdf = &gsl_simplex_log2posterior_both;

    gsl_vector_memcpy(x, minimizer->x);

    gsl_multimin_fdfminimizer_set(fdf_minimizer, &fdf_minimizer_function, x, 
                                  step_size, mode_tolerance);





    for (size_t i = 0; i != max_iterations; ++i)
//     while (gsl_multimin_test_size(simplex_size, mode_tolerance) == GSL_CONTINUE
//            && multimin_iteration < max_iterations)
    {
        gsl_multimin_fdfminimizer_iterate(fdf_minimizer);
        double f = gsl_multimin_fdfminimizer_minimum(fdf_minimizer);
        if (i % 100 == 0)
        {
            printf("function value: %20.20g, iteration: %i\n",
                   f, static_cast<int>(i));
        }
    }

    fprintf(stdout, "mode after fine method: %g\t%g\t%g\n",
            gsl_vector_get(fdf_minimizer->x, 0),
            gsl_vector_get(fdf_minimizer->x, 1),
            gsl_vector_get(fdf_minimizer->x, 2));

    for (size_t d = 0; d != ndim; ++d)
    {
        mode_point[d] = gsl_vector_get(fdf_minimizer->x, d);
    }

    //deallocate minimizer and x vector
    gsl_multimin_fminimizer_free(minimizer);
    gsl_multimin_fdfminimizer_free(fdf_minimizer);
    gsl_vector_free(x);
    gsl_vector_free(ss);
    gsl_vector_free(x_cur);
    gsl_vector_free(x_prev);
    gsl_vector_free(x_jump);

    return multimin_iteration;
}
*/


        do
        {

            metropolis.set_current_point(posterior.mode_point);

            metropolis.sample(tuning_num_points, 0, 1,
                              &proposal_mean, &proposal_variance);

            best_autocor_offset =
                best_autocorrelation_offset(metropolis.get_samples(),
                                            num_dimensions, tuning_num_points,
                                            autocor_max_offset, autocor_max_value);

            fprintf(stdout, 
                    "alpha0: %g, best_offset: %Zu, "
                    "position: %i, proposal mean: %g, "
                    "proposal variance: %g\n", 
                    dirichlet.get_alpha0(), best_autocor_offset,
                    summary._position, proposal_mean, proposal_variance);

            dirichlet.set_alpha0(dirichlet.get_alpha0() * dirichlet_spread_factor);
            dirichlet.set_alphas_from_mode(posterior.mode_point);

        }
        while (iteration < max_tuning_iterations &&
               best_autocor_offset > target_autocor_offset);



void ErrorEstimate::Log2PosteriorGradient(double const* sample_composition,
                                          double * gradient) const
{

    size_t ndim = 3;

    if (! (normalized(sample_composition, 4, 1e-10) &&
           all_positive(sample_composition, 4)))
    {
        fprintf(stderr, "Log2PosteriorGradient: invalid input.\n");
        exit(2);
    }
                
    std::fill(gradient, gradient + ndim, 0.0);

    double log2 = log(2);

    BASE_COMP_COUNTS::const_iterator bc;
    for (bc = this->base_comp_counts.begin(); 
         bc != this->base_comp_counts.end(); ++bc)
    {
        BaseComp const& evidence_based_comp = (*bc).first;
        int count = (*bc).second;
        double factor = 
            static_cast<double>(count)
            / (this->SingleObservation(evidence_based_comp.data, sample_composition)
               * log2);
        
        for (size_t d = 0; d != ndim; ++d)
        {
            gradient[d] += 
                factor
                * this->SingleObservationGradient(evidence_based_comp.data, d);
        }
    }
}



void levmar_log2posterior(double * x, double * y, int /* ndim */, 
                          int /* n_output_dim */, void * constant_params)
{
    ErrorEstimate * ee = static_cast<ErrorEstimate *>(constant_params);
    BaseComp composition(x[0], x[1], x[2], x[3]);
    if (! (normalized(composition.data, 4, 1e-10) &&
           all_positive(composition.data, 4)))
    {
        fprintf(stderr, "levmar_log2posterior: out of bounds input\n");
        exit(1);
    }

    size_t component = 0;
    BASE_COMP_COUNTS::const_iterator bc;
    for (bc = ee->base_comp_counts.begin(); 
         bc != ee->base_comp_counts.end(); ++bc)
    {
        BaseComp const& evidence_based_comp = (*bc).first;
        int count = (*bc).second;
        
        y[component] =
            log2f(ee->SingleObservation(evidence_based_comp.data, composition.data))
            * static_cast<REAL>(count) / 2.0;
        ++component;
    }
}


void gsl_simplex_log2posterior_gradient(gsl_vector const* x, void * params, 
                                        gsl_vector * gradient)
{
    GSL_SIMPLEX_PARAMS * ee_and_mode =
        static_cast<GSL_SIMPLEX_PARAMS * >(params);

    ErrorEstimate * ee = ee_and_mode->first;

    double xx[] = { 
        gsl_vector_get(x, 0), 
        gsl_vector_get(x, 1),
        gsl_vector_get(x, 2)
    };

    BaseComp composition(xx[0], xx[1], xx[2]);
    double sum = xx[0] + xx[1] + xx[2];

    if (sum > 1.0)
    {
        if (all_positive(composition.data, 4))
        {
            fprintf(stdout, "setting repel gradient.\n");
            gsl_vector_set_all(gradient, DBL_MAX);
        }
        else
        {
            for (size_t d = 0; d != 3; ++d)
            {
                double grad;
                if (xx[d] < 0){
                    grad = -DBL_MAX;
                }
                else if (xx[d] > 1)
                {
                    grad = DBL_MAX;
                }
                else
                {
                    grad = 0.0;
                }
                fprintf(stdout, "setting bounds gradient.\n");
                gsl_vector_set(gradient, d, grad);
            }
        }
    }
    else
    {
        double gtmp[3];
        ee->Log2PosteriorGradient(composition.data, gtmp);
//         fprintf(stdout, "gradient: %g\t%g\t%g\n", gtmp[0], gtmp[1], gtmp[2]);
        gsl_vector_set(gradient, 0, - gtmp[0]);
        gsl_vector_set(gradient, 1, - gtmp[1]);
        gsl_vector_set(gradient, 2, - gtmp[2]);
    }
}


void gsl_simplex_log2posterior_both(gsl_vector const* x, 
                                    void * params, 
                                    double * value, 
                                    gsl_vector * gradient)
{
    *value = gsl_simplex_log2posterior(x, params);
    gsl_simplex_log2posterior_gradient(x, params, gradient);
}




//returns P(I|b), I = image data, b = founder base
//I is currently represented by the group of all image data
//that would yield the same base call and quality score
double ErrorEstimate::DataProbability(BASE_QUAL const& image_datum,
                                      char founder_base) const
{
    double major_minor[4];
    int called_base_index = 
        ErrorEstimate::base_to_index[static_cast<size_t>(image_datum.first)];

    int founder_base_index = 
        ErrorEstimate::base_to_index[static_cast<size_t>(founder_base)];

    int quality = image_datum.second;

    double correctness_prob = 1.0 - QualityToErrorProb(quality);
    SetMajorMinorComposition(correctness_prob, called_base_index, major_minor, 4);
    return 
        this->DataPrior(image_datum) 
        * major_minor[founder_base_index]
        / base_prior[founder_base_index];
}



//Prior for the observed data
REAL ErrorEstimate::Log2DataPrior(BASE_QUAL const& image_datum) const
{
    return log2f(this->DataPrior(image_datum));
}


//computes P(I), I the image data represented by the base, qual
//pair.
double ErrorEstimate::DataPrior(BASE_QUAL const& image_datum) const
{
    int majority_bi = 
        ErrorEstimate::base_to_index[static_cast<size_t>(image_datum.first)];
    return this->data_prior[majority_bi][image_datum.second];
}



        if (false)
        {
            //create a multinomial sampling

            //create sample points
            double * multinomial_sample_points = new double[final_num_points * full_ndim];
            unsigned int * multinomial_counts = new unsigned int[full_ndim];

            sprintf(line_label, "%s\t%s\t%i\t%c\t%i\t%Zu", "MN", summary._reference, 
                    summary._position, summary._reference_base, summary._read_depth,
                    effective_depth);
        
            std::fill(multinomial_sample_points,
                      multinomial_sample_points + (final_num_points * full_ndim), 0.0);

            gsl_rng * rand_seed = gsl_rng_alloc(gsl_rng_taus);

            for (size_t s = 0; s != final_num_points; ++s)
            {
                BASE_COMP_COUNTS::const_iterator bc;
                for (bc = base_comp_counts.begin();
                     bc != base_comp_counts.end(); ++bc)
                {
                    BaseComp const& comp = (*bc).first;
                    int count = (*bc).second;
                    gsl_ran_multinomial(rand_seed, full_ndim, count, comp.data, 
                                        multinomial_counts);

                    for (size_t d = 0; d != full_ndim; ++d)
                    {
                        multinomial_sample_points[s * full_ndim + d] +=
                            static_cast<double>(multinomial_counts[d]);
                    }
                }
                //in-place normalization
                normalize(multinomial_sample_points + (s * full_ndim), full_ndim,
                          multinomial_sample_points + (s * full_ndim));

            }

            gsl_rng_free(rand_seed);

            std::vector<double *> sample_points_sortable(final_num_points);
        
            for (size_t i = 0; i != final_num_points; ++i)
            {
                sample_points_sortable[i] = 
                    multinomial_sample_points + (i * full_ndim);
            }

            double fake_mode[4] = { -1, -1, -1, -1 };

            print_quantiles(posterior_output_fh, & sample_points_sortable, 
                            fake_mode,
                            line_label, dimension_labels, "+", quantiles,
                            num_quantiles, full_ndim);

            delete multinomial_sample_points;
            delete multinomial_counts;
        }




        int m;
        lbfgsstate state;
        lbfgsreport rep;
        ap::real_1d_array comp_array;
        double current_comp[truncated_ndim];
        double current_grad[truncated_ndim];
        double adjusted_comp[truncated_ndim];
        double adjusted_grad[truncated_ndim];


        //
        // Function minimized:
        //     F = exp(x-1) + exp(1-x) + (y-x)^2
        // N = 2 - task dimension
        // M = 1 - build tank-1 model
        //
        m = 1;
        comp_array.setlength(truncated_ndim);
        comp_array(0) = 0.25;
        comp_array(1) = 0.25;
        comp_array(2) = 0.25;
//         comp_array(3) = 0.25;
        double eps_gradient = 1e-10;
        double eps_funcval = 1e-10;
        double eps_coord = 1e-10;

        size_t optim_ndim = truncated_ndim;

        minlbfgs(optim_ndim, m, comp_array, eps_gradient, 
                 eps_funcval, eps_coord, 0, 0, state);

        //while(minlbfgsiteration(state))
        while(1)
        {
            minlbfgsiteration(state);
            
            for (size_t d = 0; d != optim_ndim; ++d)
            {
                current_comp[d] = state.x(d);
            }
            current_comp[3] = 1.0 - current_comp[0] - current_comp[1] - current_comp[2];

            if (all_positive(current_comp, full_ndim) &&
                normalized(current_comp, full_ndim, 1e-10))
            {
//                 normalize(current_comp, full_ndim, current_comp);

                error_estimate.log2_posterior_gradient(current_comp, current_grad);
                for (size_t d = 0; d != optim_ndim; ++d)
                {
                    state.g(d) = - current_grad[d];
                }
                state.f = - error_estimate.Log2Posterior(current_comp);
            }
            else
            {
                
                //distance
                double dist =
                    (current_comp[0] - 0.25) * (current_comp[0] - 0.25)
                    + (current_comp[1] - 0.25) * (current_comp[1] - 0.25)
                    + (current_comp[2] - 0.25) * (current_comp[2] - 0.25);

                state.f = expf(dist) + 1e20;
                state.g(0) = expf(dist) * 2.0 * (current_comp[0] - 0.25);
                state.g(1) = expf(dist) * 2.0 * (current_comp[1] - 0.25);
                state.g(2) = expf(dist) * 2.0 * (current_comp[2] - 0.25);

            }

        }
        minlbfgsresults(state, comp_array, rep);





        //
        // output results
        //
        printf("F = exp(x-1) + exp(1-x) + (y-x)^2\n");
//         printf("X = %4.2lf (should be 1.00)\n", static_cast<double>(comp_array(0)));
        printf("mode point: (%10.8g\t%10.8g\t%10.8g\n",
               comp_array(0), comp_array(1), comp_array(2));


//         SmartPtr<TNLP> nlp_posterior = new IpoptPosterior(error_estimate, 3);
        
//         Ipopt::SmartPtr<Ipopt::IpoptApplication> posterior_ipopt_app = 
//             new Ipopt::IpoptApplication();
        
        // Change some options
        // Note: The following choices are only examples, they might not be
        //       suitable for your optimization problem.
//         posterior_ipopt_app->Options()->SetNumericValue("tol", 1e-7);
//         posterior_ipopt_app->Options()->SetStringValue("mu_strategy", "adaptive");
//         posterior_ipopt_app->Options()->SetStringValue("output_file", "ipopt.out");
//         posterior_ipopt_app->Options()->SetStringValue("hessian_approximation", "limited-memory");
//         posterior_ipopt_app->Options()->SetStringValue("linear_solver", "mumps");
//         posterior_ipopt_app->Options()->SetStringValue("derivative_test", "first-order");
        
        // The following overwrites the default name (ipopt.opt) of the
        // options file
        // posterior_ipopt_app->Options()->SetStringValue("option_file_name", "hs071.opt");

        // Intialize the IpoptApplication and process the options
//         Ipopt::ApplicationReturnStatus status;
        
//         status = posterior_ipopt_app->Initialize();

//         if (status != Ipopt::Solve_Succeeded) {
//             printf("\n\n*** Error during initialization!\n");
//             return (int) status;
//         }

//         // Ask Ipopt to solve the problem
//         status = posterior_ipopt_app->OptimizeTNLP(nlp_posterior);

//         if (status == Ipopt::Solve_Succeeded) {
//             printf("\n\n*** The problem solved!\n");
//         }
//         else {
//             printf("\n\n*** The problem FAILED!\n");
//         }


    gsl_multimin_fdfminimizer * sec_minimizer =
        gsl_multimin_fdfminimizer_alloc(gsl_multimin_fdfminimizer_vector_bfgs, 
                                        sphere_ndim);

    gsl_multimin_fdfminimizer_set(sec_minimizer, &fdf_minimizer_function, first_minimizer->x, 1e-10, 0.1);

    double current_neg_mode = Transformation::log2_neg_posterior(first_minimizer->x, & passing_params);

    passing_params.current_mode = - current_neg_mode;

    multimin_iteration = 0;
    while (gsl_multimin_test_gradient(gsl_multimin_fdfminimizer_gradient(sec_minimizer), 
                                      max_gradient_norm) == GSL_CONTINUE &&
           multimin_iteration < max_iterations)
    {
        gsl_vector * gradient = gsl_multimin_fdfminimizer_gradient(sec_minimizer);
        Transformation::probability_transform(sec_minimizer->x, mode_point, mode_gradient);

        log2_neg_gradient_numerical(sec_minimizer->x, & passing_params,
                                    numerical_gradient_epsilon, numerical_gradient);

        if (multimin_iteration % 1 == 0)
        {
            printf(
                   "gslfdf2: " "%20.18g\t" 
                   "%20.18g\t" "%20.18g\t" "%20.18g\t"
                   "%20.18g\t" "%20.18g\t" "%20.18g\t"
                   "%20.18g\t" "%20.18g\t" "%20.18g\t" "%20.18g\t"
                   "%20.18g\t" "%20.18g\t" "%20.18g\n"
                   ,
                   Transformation::log2_neg_posterior(sec_minimizer->x, & passing_params),
                   gsl_vector_get(gradient, 0),
                   gsl_vector_get(numerical_gradient, 0),
                   gsl_vector_get(gradient, 1),
                   gsl_vector_get(numerical_gradient, 1),
                   gsl_vector_get(gradient, 2),
                   gsl_vector_get(numerical_gradient, 2),
                   mode_point[0], mode_point[1], mode_point[2], mode_point[3],
                   gsl_vector_get(sec_minimizer->x, 0),
                   gsl_vector_get(sec_minimizer->x, 1),
                   gsl_vector_get(sec_minimizer->x, 2)
                   );
        }

        gsl_multimin_fdfminimizer_iterate(sec_minimizer);

//         double norm = gsl_blas_dnrm2(gradient);

        ++multimin_iteration;
    }

    gsl_vector_memcpy(x, sec_minimizer->x);


    if (false)
    {
        //print out a bunch of grid points centered around the mode
        double total_width = 1e-2;
        size_t num_grid_points = 5;
        double step = total_width / static_cast<double>(num_grid_points);
        double half_width = total_width / 2.0;
        
        double func_val;
        double comp[4];
        double comp_gradient[4][3];
        
        gsl_vector * sphere_gradient = gsl_vector_alloc(sphere_ndim);
        
        //explore the space, setting each of the total_widths by line search
        gsl_vector_set(x, 0, gsl_vector_get(first_minimizer->x, 0) - half_width);
        gsl_vector_set(x, 1, gsl_vector_get(first_minimizer->x, 1) - half_width);
        gsl_vector_set(x, 2, gsl_vector_get(first_minimizer->x, 2) - half_width);
        
        for (size_t d0 = 0; d0 != num_grid_points; ++d0)
        {
            gsl_vector_set(x, 0, gsl_vector_get(first_minimizer->x, 0) - (total_width / 2.0) + (d0 * step));
            for (size_t d1 = 0; d1 != num_grid_points; ++d1)
            {
                gsl_vector_set(x, 1, gsl_vector_get(first_minimizer->x, 1) - (total_width / 2.0) + (d1 * step));
                for (size_t d2 = 0; d2 != num_grid_points; ++d2)
                {
                    gsl_vector_set(x, 2, gsl_vector_get(first_minimizer->x, 2) - (total_width / 2.0) + (d2 * step));

                    Transformation::log2_neg_posterior_and_gradient(x, & passing_params,
                                                                    & func_val, sphere_gradient);

                    Transformation::probability_transform(x, comp, comp_gradient);

                    printf(
                           "grid: " "%20.10f\t" 
                           "%20.10f\t" "%20.10f\t" "%20.10f\t"
                           "%20.10f\t" "%20.10f\t" "%20.10f\t" "%20.10f\t"
                           "%20.10f\t" "%20.10f\t" "%20.10f\n"
                           ,
                           func_val,
                           gsl_vector_get(sphere_gradient, 0),
                           gsl_vector_get(sphere_gradient, 1),
                           gsl_vector_get(sphere_gradient, 2),
                           comp[0], comp[1], comp[2], comp[3],
                       
                           gsl_vector_get(x, 0),
                           gsl_vector_get(x, 1),
                           gsl_vector_get(x, 2)
                           );
                }
                printf("\n");
            }
            printf("\n");
        }
    }






/*
  Probability of observing a set of reported nucleotide probabilities
  given a particular sample composition.
  Although the reported nucleotide probabilities are themselves probabilities,
  they are treated here as observations.  They are also used as probabilities.
*/
REAL ErrorEstimate::SingleObservation(double const* reported_nuc_probs,
                                      double const* sample_comp) const
{
    return
        (reported_nuc_probs[0] * sample_comp[0] +
         reported_nuc_probs[1] * sample_comp[1] +
         reported_nuc_probs[2] * sample_comp[2] +
         reported_nuc_probs[3] * sample_comp[3])
        * 1.0;
    //we are not keeping track of the data prior here since it remains constant
    //over the sample composition

}


double ErrorEstimate::SingleObservationGradient(double const* reported_nuc_probs,
                                                size_t gradient_dimension) const
{
    return reported_nuc_probs[gradient_dimension];
}


BASE_QUAL_COUNTS TallyBaseQualCounts(char const* bases, 
                                     char const* quality_codes,
                                     int num_bases_quals)
{

    BASE_QUAL_COUNTS base_qual_counts;

    int qual;
    for (int q = 0; q < num_bases_quals; ++q)
    {
        
        if (bases[q] == '*')
        {
            //not a real base.
            continue;
        }

        qual = QualityCodeToQuality(quality_codes[q]);
        assert(qual <= 34);

        std::pair<char, int> base_qual(bases[q], qual);

        if (base_qual_counts.find(base_qual) == base_qual_counts.end())
        {
            base_qual_counts.insert(std::make_pair(base_qual, 0));
        }
        base_qual_counts[base_qual]++;
    }
    return base_qual_counts;
}


//assume x is somewhere in the hypercube
double gsl_simplex_hypercube_posterior(gsl_vector const* x, void * params)
{

    PosteriorParams * pp = static_cast<PosteriorParams *>(params);

    double xx[] = { 
        gsl_vector_get(x, 0), 
        gsl_vector_get(x, 1),
        gsl_vector_get(x, 2) 
    };
    double xc[3];

    if (! within_hypercube(xx))
    {
        return DBL_MAX;
    }

    pp->ee->contract_from_hypercube(xx, xc);

    BaseComp composition(xc[0], xc[1], xc[2]);
        
    return -1.0 * pp->ee->Log2Posterior(composition.data);
}


double gsl_simplex_log2posterior(gsl_vector const* x, void * params)
{

    GSL_SIMPLEX_PARAMS * ee_and_mode =
        static_cast<GSL_SIMPLEX_PARAMS * >(params);

    ErrorEstimate * ee = ee_and_mode->first;

    double xx[] = { 
        gsl_vector_get(x, 0), 
        gsl_vector_get(x, 1),
        gsl_vector_get(x, 2),
        0.0
    };

    xx[3] = 1.0 - xx[0] - xx[1] - xx[2];
    double val;

    if (! (normalized(xx, 4, 1e-10) &&
           all_positive(xx, 4)))
    {
        val = DBL_MAX;
    }
    else
    {
        val = -1.0 * ee->Log2Posterior(xx);
    }

//     fprintf(stdout, "point: %g\t%g\t%g\tval: %g\n", xx[0], xx[1], xx[2], val);
    return val;
}



//Count the number of base_of_interest in the bases
template<typename T>
int CountElement(T const* elements, int num_elements, T element_of_interest)
{
    int num_elements_of_interest = 0;
    for (int b = 0; b < num_elements; ++b)
    {
        num_elements_of_interest +=
            (elements[b] == element_of_interest) ? 1 : 0;
    }
    return num_elements_of_interest;
}


//Return a pair of <majority_base, count> for the bases
template<typename T>
std::pair<char, int> MajorityElement(T const* elements, int size, 
                                     T const* alphabet, int alphabet_size)
{

    //determine majority base
    int * counts = new int[alphabet_size];

    int max_count_index = 0;
    int max_count_value = 0;

    for (int b = 0; b < alphabet_size; ++b)
    {
        counts[b] = CountElement(elements, size, alphabet[b]);
        if (counts[b] > max_count_value)
        {
            max_count_value = counts[b];
            max_count_index = b;
        }
    }
    T mb = elements[max_count_index];
    delete counts;

    return std::make_pair(mb, max_count_value);
}

typedef std::pair<REAL, size_t> COORDINATE;

struct less_coordinate
{
    bool operator()(COORDINATE const& c1,
                    COORDINATE const& c2) const;
};





//calculate the mean and covariance from the set of weighted samples
//storing them in mu and sigma
//covariance is a single array 
void CalculateMeanCovariance(WEIGHTED_SAMPLE_MAP const& samples, 
                             size_t const ndim,
                             double * mu, double * covariance)
{
    WEIGHTED_SAMPLE_MAP::const_iterator sample_iter;

    std::vector<double> coord;

    //collapse sample coords into flat 'coord'
    for (sample_iter = samples.begin(); 
         sample_iter != samples.end(); ++sample_iter)
    {
        WeightedSample const& ws = (*sample_iter).second;
        for (size_t w = 0; w != static_cast<size_t>(ws.weight); ++w)
        {
            for (size_t d = 0; d != ndim; ++d)
            {
                coord.push_back(ws.x[d]);
            }
        }
    }

    //the number of samples, counting the total weight
    size_t num_effective_samples = coord.size() / ndim;

    //update the mu's
    for (size_t d = 0; d != ndim; ++d)
    {
        mu[d] = gsl_stats_mean(& (*coord.begin()) + d, ndim, 
                               num_effective_samples);
    }
    
    //update the covariances
    for (size_t d1 = 0; d1 != ndim; ++d1)
    {
        for (size_t d2 = 0; d2 != ndim; ++d2)
        {
            covariance[d1 * ndim + d2] = 
                gsl_stats_covariance_m(&coord[d1], ndim,
                                       &coord[d2], ndim,
                                       num_effective_samples,
                                       mu[d1], mu[d2]);
        }
    }
}


//initialize the relevant cdf[] field in the samples,
//and leave <samples> sorted according to <keep_dimension>
void MarginalCumulativeDistribution(size_t keep_dimension,
                                    std::vector<WeightedSample *> * samples)
{

    std::sort((*samples).begin(), (*samples).end(), SortCoordinate(keep_dimension));
    
    REAL cumulant = 0.0;
    for (size_t s = 0; s != (*samples).size(); ++s)
    {
        cumulant += (*samples)[s]->weight;
        (*samples)[s]->cdf[keep_dimension] = cumulant;
    }

}



size_t ErrorEstimate::find_mode_point2(double max_gradient_norm, 
                                       size_t max_iterations,
                                       double * mode_point) const
{
    lbfgsstate state;
//     lbfgsreport rep;

    size_t sphere_ndim = 3;

    int m = 1;

    double eps_gradient = 1e-10;
    double eps_funcval = 1e-10;
    double eps_coord = 1e-10;

    gsl_vector * sphere_coords = gsl_vector_alloc(sphere_ndim);
    gsl_vector_set(sphere_coords, 0, 0.25);
    gsl_vector_set(sphere_coords, 1, 0.25);
    gsl_vector_set(sphere_coords, 2, 0.25);

    ap::real_1d_array start_point;
    start_point.setlength(sphere_ndim);
    start_point(0) = 0.25;
    start_point(1) = 0.25;
    start_point(2) = 0.25;

    gsl_vector * sphere_gradient = gsl_vector_alloc(sphere_ndim);

    minlbfgs(sphere_ndim, m, start_point, eps_gradient, 
             eps_funcval, eps_coord, 0, 0, state);


    Transformation::PassingParams passing_params = { this, 0.0 };

    double comp_gradient[4][3];


    //while(minlbfgsiteration(state))
    //while(1)
    for (size_t iter = 0; iter != max_iterations; ++iter)
    {
        minlbfgsiteration(state);

        for (size_t d = 0; d != sphere_ndim; ++d)
        {
            gsl_vector_set(sphere_coords, d, state.x(d));
        }
        Transformation::log2_neg_posterior_and_gradient(sphere_coords, & passing_params,
                                                        & state.f, sphere_gradient);
                                                            
        for (size_t d = 0; d != sphere_ndim; ++d)
        {
            state.g(d) = gsl_vector_get(sphere_gradient, d);
        }

        Transformation::probability_transform(sphere_coords, mode_point, comp_gradient);
        printf(
               "alglib: "
               "%20.18g\t" "%20.18g\t" "%20.18g\t" "%20.18g\t"
               "%20.18g\t" "%20.18g\t" "%20.18g\t" "%20.18g\t"
               "%20.18g\t" "%20.18g\t" "%20.18g\n"
               ,
               Transformation::log2_neg_posterior(sphere_coords, & passing_params),
               gsl_vector_get(sphere_gradient, 0),
               gsl_vector_get(sphere_gradient, 1),
               gsl_vector_get(sphere_gradient, 2),
               mode_point[0], mode_point[1], mode_point[2], mode_point[3],
               gsl_vector_get(sphere_coords, 0),
               gsl_vector_get(sphere_coords, 1),
               gsl_vector_get(sphere_coords, 2)
               );

    }

//     ap::real_1d_array final_sphere_coords;
//     final_sphere_coords.set_length(sphere_ndim);
//     minlbfgsresults(state, final_sphere_coords, rep);

    Transformation::probability_transform(sphere_coords, mode_point, comp_gradient);

    gsl_vector_free(sphere_gradient);
    return 100;
}


//Group identical numbers of base compositions
OBSERVED_DATA_COUNTS TallyDataCounts(std::vector<ObservedData const*> const& compositions)
{
    OBSERVED_DATA_COUNTS counts;
    std::vector<ObservedData const*>::const_iterator cit;
    for (cit = compositions.begin(); cit != compositions.end(); ++cit)
    {
        if (counts.find((*cit)) == counts.end())
        {
            counts.insert(std::make_pair((*cit), 0));
        }
        counts[(*cit)]++;
    }
    return counts;
}


//Group identical numbers of base compositions with strand
void TallyDataCountsStrand(std::vector<ObservedData const*> const& compositions,
                           OBSERVED_DATA_COUNTS * pos_strand_counts,
                           OBSERVED_DATA_COUNTS * neg_strand_counts)
{

    for (size_t read = 0; read != compositions.size(); ++read)
    {
        ObservedData const* comp = compositions[read];

        OBSERVED_DATA_COUNTS & chosen_counts = comp.strand == POS_STRAND 
            ? *pos_strand_counts : *neg_strand_counts;

        if (chosen_counts.find(comp) == chosen_counts.end())
        {
            chosen_counts.insert(std::make_pair(comp, 0));
        }
        chosen_counts[comp]++;
    }
}




//     BaseQualStrand(REAL _a, REAL _c, REAL _g, REAL _t, DNAStrand _strand) : 
//         A(_a), C(_c), G(_g), T(_t), strand(_strand)
//     { 
//         data[0] = A;
//         data[1] = C;
//         data[2] = G;
//         data[3] = T;
//         majority_element =
//             std::distance(data, std::max_element(data, data + 4));
//         quality_score = error_prob_to_quality(1.0 - data[majority_element]);
//     }

//     BaseQualStrand(REAL _a, REAL _c, REAL _g, DNAStrand _strand) : 
//         A(_a), C(_c), G(_g), T(1.0 - _a - _c - _g), strand(_strand)
//     { 
//         data[0] = A;
//         data[1] = C;
//         data[2] = G;
//         data[3] = T;
//         majority_element =
//             std::distance(data, std::max_element(data, data + 4));
//         quality_score = error_prob_to_quality(1.0 - data[majority_element]);
//     }


    /*     { */
/*         double error_prob = QualityToErrorProb(quality_score); */
/*         double other = error_prob / 3.0; */
/*         std::fill(data, data+4, other); */
/*         data[majority_element] = 1.0 - error_prob; */
/*         A = data[0]; */
/*         C = data[1]; */
/*         G = data[2]; */
/*         T = data[3]; */
/*     } */
    
//     bool valid(REAL epsilon) const
//     {
//         return positive(epsilon) && fabsl(norm() - 1.0) <= epsilon;
//     }
    
//     BaseQualStrand() : A(0.25), C(0.25), G(0.25), T(0.25), strand(POS_STRAND)
//     { 
//         data[0] = A;
//         data[1] = C;
//         data[2] = G;
//         data[3] = T;
//         majority_element =
//             std::distance(data, std::max_element(data, data + 4));
//         quality_score = error_prob_to_quality(1.0 - data[majority_element]);
//     }
    
    
//     bool positive(REAL epsilon) const
//     {
//         return A >= -epsilon
//             && C >= -epsilon
//             && G >= -epsilon
//             && T >= -epsilon;
//     }

//     REAL norm() const
//     {
//         return A + C + G + T;
//     }


//initializes the data JPD from a generic jpd
void ErrorEstimate::Initialize(OBSERVED_DATA_JPD const& data_jpd)
{

    OBSERVED_DATA_JPD::const_iterator data_iter;
    size_t count = 0;
    double norm[NBASES];
    double norm_total;

    std::fill(norm, norm + NBASES, 0.0);

    for (data_iter = data_jpd.begin(); data_iter != data_jpd.end(); ++data_iter)
    {
        size_t & ct[4] = (*data_iter).second;
        this->observed_data_order[(*data_iter).first] = count++;
        
        for (size_t base_index = 0; base_index != NBASES; ++base_index)
        {
            norm[base_index] += ct[base_index];
        }
    }

    norm_total = std::accumulate(norm, norm + NBASES, 0.0);

    size_t num_distinct_data = data_jpd.size();

    for (size_t base_index = 0; base_index != NBASES; ++base_index)
    {
        this->observed_data_jpd[base_index] = new double[num_distinct_data];
        for (data_iter = data_jpd.begin(); data_iter != data_jpd.end(); ++data_iter)
        {
            size_t data_index = this->observed_data_order[*data_iter];
            observed_data_jpd[base_index][data_index] = 
                static_cast<double>((*data_iter).second[base_index]) / norm_total;
        }
    }
    
    for (size_t base_index = 0; base_index != NBASES; ++base_index)
    {
        this->founder_base_prob[base_index] =
            static_castnorm[base_index] / norm_total;
    }
    
}


struct less_base_comp
{
    bool operator()(BaseComp const& a,
                    BaseComp const& b) const
    {
        return 
            a.A < b.A
            || (a.A == b.A 
                && (a.C < b.C 
                    || (a.C == b.C 
                        && (a.G < b.G
                            || (a.G == b.G
                                && (a.T < b.T
                                    || (a.T == b.T
                                        && a.strand < b.strand)))))));
    }
};


typedef std::map<BaseComp, int, less_base_comp> BASE_COMP_COUNTS;


/*

  Illumina example: the tuple (calle_base, quality, strand) 

  a generic class for observable data.  describes all observed data
  for one base position of one read.

   JPD P(founder_base, observed_data).

   also, provides an index for quick lookup

*/



//expand each raw count according to the definition of Phred-quality and 
//assumption of uniformity over the three error bases.
void set_data_from_raw_counts(std::map<std::string, size_t> const& name_map,
                              double const* raw_counts,
                              std::set<BaseQualStrandReader::Datum> const& raw_counts,
                              NucleotideStats * stats)
{
    std::set<BaseQualStrandReader::Datum>::const_iterator data_iter;
    double total_frequency = 0.0;

    for (name_iter = name_map.begin(); name_iter != name_map.end(); ++name_iter)
    for (data_iter = raw_counts.begin(); data_iter != raw_counts.end(); ++data_iter)
    {
        
        BaseQualStrandReader::Datum datum =
            get_datum_from_name((*name_iter).first);

        size_t index = (*name_iter).second;
        double frequency = raw_counts[index];

        BaseQualStrandReader::Datum const& datum = (*data_iter);
        total_frequency += datum.frequency;

        double dist[4];
        double error_prob = QualityToErrorProb(datum.quality);
        double other_prob = error_prob / 3.0;

        size_t basecall_index = 
            Nucleotide::base_to_index[static_cast<int>(datum.called_base)];

        //in this application we forbid 'N'.  only ACGT
        assert(basecall_index < 4);
        
        std::fill(dist, dist+4, other_prob * datum.frequency);
        dist[basecall_index] = datum.frequency * (1.0 - error_prob);

        for (size_t b = 0; b != 4; ++b)
        {
            stats->complete_jpd[b][datum.order_index] = dist[b];
        }

        stats->raw_counts[datum.order_index] = datum.frequency;
        stats->name_mapping[datum.name()] = datum.order_index;
        stats->index_mapping[datum.order_index] = datum.name();
    }

    size_t D = stats->num_distinct_data;
    for (size_t b = 0; b != 4; ++b)
    {
        for (size_t di = 0; di != D; ++di)
        {
            stats->complete_jpd[b][di] /= total_frequency;
        }
        stats->founder_base_marginal[b] =
            std::accumulate(stats->complete_jpd[b],
                            stats->complete_jpd[b] + D, 0.0);
    }

    for (size_t b = 0; b != 4; ++b)
    {
        for (size_t di = 0; di != D; ++di)
        {
            stats->founder_base_likelihood[b][di] =
                stats->complete_jpd[b][di]
                / stats->founder_base_marginal[b];
        }
    }

}


REAL epsilon_gradient(double * x, double * mode_point, double a, double c)
{
    double * m = mode_point;
    double dist_sq = (x[0]-m[0]) * (x[0]-m[0]) + (x[1]-m[1]) * (x[1]-m[1]) + (x[2]-m[2]) * (x[2]-m[2]); 
    return a / (sqrt(dist_sq) + c);
}


SliceSampleComposition : SliceSampleComposition.o error_estimate.o dirichlet.o pileup_tools.o hilbert.o \
	tools.o slice_sampling.o sampling.o integrands.o
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LDFLAGS)



//sets the output distribution to uniform for all non-major indices,
//and m for the i.e {0.6, 0.1, 0.1, 0.1, 0.1}
void SetMajorMinorComposition(double major_probability, int major_index,
                              double * output_distribution, int size)
{
    for (int i = 0; i != size; ++i)
    {
        output_distribution[i] = (i == major_index) ?
            major_probability :
            (1.0 - major_probability) / (size - 1);
    }
}


//simulate the process of measuring an actual base
std::pair<char, int> SimulateBaseMeasure(gsl_rng * rand_gen,
                                         int const* quality_counts, 
                                         int total_counts, char base)
{
    int quality = 
        SampleDiscreteDistribution(rand_gen, quality_counts, total_counts);

    double good_read_prob = 1.0 - QualityToErrorProb(quality);
    
    int major_index = Nucleotide::base_to_index[static_cast<int>(base)];

    double miscall_distribution[4];
    SetMajorMinorComposition(good_read_prob, major_index, miscall_distribution, 4);

    int measured_base_index = 
        SampleDiscreteDistribution(rand_gen, miscall_distribution, 1.0);

    return std::make_pair(Nucleotide::bases_upper[measured_base_index],
                          quality);
}


//         double log_corner_values[full_ndim];
//         for (size_t d = 0; d != full_ndim; ++d)
//         {
//             double corner_point[full_ndim];
//             std::fill(corner_point, corner_point + full_ndim, 1e-10);
//             corner_point[d] = 1.0 - (3.0 * 1e-10);
//             log_corner_values[d] = posterior.log_pdf(corner_point);
//         }


        if (0)
            //if (1)
        {
            gsl_vector * g = gsl_vector_alloc(3);
            gsl_vector * r = gsl_vector_alloc(3);
            double value;
            double xlow[3];
            double comp[4];
            double grid_spacing = 1;
            double npoints = 5;
            Transformation::SigmoidVals sigmoid_vals[3];

            Transformation::PassingParams pass_params = { & model, 0.0 };
            for (int dim = 0; dim != 3; ++dim)
            {
                xlow[dim] = - ((npoints - 1) / 2.0) * grid_spacing;
            }

            for (size_t x = 0; x != npoints; ++x)
            {
                double xc = xlow[0] + x * grid_spacing;
                gsl_vector_set(r, 0, xc);
                for (size_t y = 0; y != npoints; ++y)
                {
                    double yc = xlow[1] + y * grid_spacing;
                    gsl_vector_set(r, 1, yc);
                    for (size_t z = 0; z != npoints; ++z)
                    {
                        double zc = xlow[2] + z * grid_spacing;
                        gsl_vector_set(r, 2, zc);

                        Transformation::log_neg_posterior_value_and_gradient(r, 
                                                                             static_cast<void *>(& pass_params),
                                                                             &value, g);
                        double x[3];
                        x[0] = gsl_vector_get(r, 0);
                        x[1] = gsl_vector_get(r, 1);
                        x[2] = gsl_vector_get(r, 2);
                        sigmoid_value_and_gradient(x, sigmoid_vals);
                        sigmoid_composition(sigmoid_vals, comp);
                        
                        fprintf(stdout, 
                                "%-10.8g\t%-10.8g\t%-10.8g\t"
                                "%-10.8g\t%-10.8g\t%-10.8g\t"
                                "%-10.8g\t"
                                "%-10.8g\t%-10.8g\t%-10.8g\t%-10.8g\n",
                                gsl_vector_get(r, 0),
                                gsl_vector_get(r, 1),
                                gsl_vector_get(r, 2),
                                - gsl_vector_get(g, 0),
                                - gsl_vector_get(g, 1),
                                - gsl_vector_get(g, 2),
                                - value,
                                comp[0], comp[1], comp[2], comp[3]);
                    }
                }
            }
            gsl_vector_free(g);
            gsl_vector_free(r);
        }
