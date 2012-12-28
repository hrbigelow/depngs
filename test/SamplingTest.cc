#include <vector>
#include <utility>


#include "henry/tools.h"
#include "henry/slice_sampling.h"
#include "henry/error_estimate.h"
#include "henry/integrands.h"
#include "henry/metropolis.h"
#include "henry/sampling.h"
#include "henry/stats_tools.h"

int main(int argc, char ** argv)
{

    char * target_mu_file = argv[1];
    char * target_covariance_file = argv[2];
    char * quantiles_file = argv[3];
    int log2_integrand = atoi(argv[4]);
    size_t number_estimation_points = static_cast<size_t>(atof(argv[5]));
    size_t number_sample_points = static_cast<size_t>(atof(argv[6]));


    printf ("target_mu_file: %s\n"
            "target_covariance_file: %s\n"
            "quantiles_file: %s\n"
            "log2_integrand: %i\n"
            "number_estimation_points: %Zu\n"
            "number_sample_points: %Zu\n"
            ,

            target_mu_file, 
            target_covariance_file, 
            quantiles_file,
            log2_integrand,
            number_estimation_points,
            number_sample_points

);



    //parse mass fractions file
    double quantiles[1000];
    size_t num_quantiles = 
        ParseNumbersFile(quantiles_file, quantiles);

    double target_mu[1000];
    size_t num_dimensions = ParseNumbersFile(target_mu_file, target_mu);

    double target_covariance[1000];
    size_t num_covariances = ParseNumbersFile(target_covariance_file, target_covariance);

    int num_bits_per_dim = 40;
    size_t initial_sampling_range = num_bits_per_dim * num_dimensions;

    if ((num_dimensions * num_dimensions) != num_covariances)
    {
        fprintf(stderr, 
                "Number of mu's (%Zu) squared not equal to number of covariances (%Zu)\n",
                num_dimensions, num_covariances);
        exit(1);
    }


    Gaussian target_density_function(target_mu, target_covariance, num_dimensions, true);

    target_density_function.Init();

    size_t range_delta = 1;

    SliceSampling slice_sampling(num_dimensions, num_bits_per_dim, 
                                 log2_integrand, range_delta);

    slice_sampling.Initialize();

    size_t autocor_max_offset = 50;
    double autocor_max_value = 0.1;

    double * estimation_points = new double[num_dimensions * number_estimation_points];
    double * sample_points = new double[num_dimensions * number_sample_points];


    slice_sampling.sample(&target_density_function, target_mu,
                          initial_sampling_range, 1,
                          estimation_points,
                          number_estimation_points);

    size_t best_autocor_offset =
        best_autocorrelation_offset(estimation_points, num_dimensions, 
                                    number_estimation_points, 
                                    autocor_max_offset, autocor_max_value);

    fprintf(stdout, "best offset: %Zu\n", best_autocor_offset);

    slice_sampling.sample(&target_density_function, target_mu,
                          initial_sampling_range, best_autocor_offset,
                          sample_points,
                          number_sample_points);


    std::vector<double *> sample_points_sortable(number_sample_points);

    for (size_t i = 0; i != number_sample_points; ++i)
    {
        sample_points_sortable[i] = sample_points + (i * num_dimensions);
    }

    print_cdf_comparison(stdout, &target_density_function,
                         & sample_points_sortable, quantiles, 
                         num_quantiles, num_dimensions);
    
    delete estimation_points;
    delete sample_points;

    return 0;

}
