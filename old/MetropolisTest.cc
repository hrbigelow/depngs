#include "henry/integrands.h"
#include "henry/error_estimate.h"


int main(int argc, char ** argv)
{

    char * target_mu_file = argv[1];
    char * target_sigma_file = argv[2];
    size_t number_sample_points = static_cast<size_t>(atof(argv[3]));
    char * quantile_points_file = argv[4];
    char * grid_points_file = argv[5];

    printf ("target_mu_file: %s\n"
            "target_sigma_file: %s\n"
            "number_sample_points: %Zu\n"
            "quantile_points_file: %s\n"
            "grid_points_file: %s\n"
            ,

            target_mu_file, 
            target_sigma_file,
            number_sample_points,
            quantile_points_file,
            grid_points_file
            );


    WEIGHTED_SAMPLE_MAP random_samples;

    double target_mu[1000];
    size_t num_dimensions = ParseNumbersFile(target_mu_file, target_mu);

    double target_sigma[1000];
    size_t num_sigmas = ParseNumbersFile(target_sigma_file, target_sigma);

    double grid_points[3];
    ParseNumbersFile(grid_points_file, grid_points);

    double min_grid_point = grid_points[0];
    double max_grid_point = grid_points[1];
    double grid_step = grid_points[2];

    double quantile_points[100000];
    size_t number_quantiles = ParseNumbersFile(quantile_points_file, quantile_points);

    if ((num_dimensions * num_dimensions) != num_sigmas)
    {
        fprintf(stderr, 
                "Number of mu's (%Zu) squared not equal to number of sigmas (%Zu)\n",
                num_dimensions, num_sigmas);
        exit(1);
    }

    Gaussian gaussian(target_mu, target_sigma, num_dimensions);

    gaussian.Init();

    double * sample_x = new double[num_dimensions];

    std::pair<WEIGHTED_SAMPLE_MAP::iterator, bool> sample_insertion;

    //test the sampling routine.
    for (size_t s = 0; s != number_sample_points; ++s)
    {
        gaussian.sample(sample_x);

        WeightedSample sample(num_dimensions, sample_x, 0.0, 0.0);
        sample_insertion = 
            random_samples.insert(std::make_pair(sample, sample));
        (*sample_insertion.first).second.weight += 1.0;
        
    }
    delete sample_x;

    double * estimated_mu = new double[num_dimensions];
    double * estimated_sigma = new double[num_dimensions * num_dimensions];

    CalculateMeanCovariance(random_samples, num_dimensions,
                            estimated_mu, estimated_sigma);

    Gaussian estimated_gaussian(estimated_mu, estimated_sigma, num_dimensions);
    gaussian.print_params(stdout);
    estimated_gaussian.print_params(stdout);

    delete estimated_mu;
    delete estimated_sigma;


    //test the operator(), cdf, and cdf_inv routines
    WEIGHTED_SAMPLE_MAP grid_samples;
    size_t num_grid_points = 
        pow((max_grid_point - min_grid_point) / grid_step, 
            num_dimensions);
    
    double * grid_point = new double[num_dimensions];
    std::fill(grid_point, grid_point + num_dimensions, min_grid_point);
    for (size_t s = 0; s != num_grid_points; ++s)
    {
        for (size_t d = 0; d != num_dimensions; ++d)
        {
            grid_point[d] += grid_step;
            if (grid_point[d] > max_grid_point)
            {
                grid_point[d] = min_grid_point;
            }
            else
            {
                break;
            }
        }
        double func_val = gaussian(grid_point, num_dimensions);
        WeightedSample sample(num_dimensions, grid_point, func_val, 
                              func_val * pow(grid_step, num_dimensions) );

        sample_insertion = 
            grid_samples.insert(std::make_pair(sample, sample));
        
    }
    delete grid_point;


    PrintCDFComparison(stdout, &gaussian, grid_samples,
                       quantile_points, number_quantiles,
                       num_dimensions);

    fprintf(stdout, "****************\n*****************\n");

    PrintCDFComparison(stdout, &gaussian, random_samples,
                       quantile_points, number_quantiles,
                       num_dimensions);


    if (false)
    {
        std::vector<WeightedSample *> grid_samples_vec;
        WEIGHTED_SAMPLE_MAP::iterator sample_iter;
        
        for (sample_iter = grid_samples.begin();
             sample_iter != grid_samples.end(); ++sample_iter)
        {
            WeightedSample * ws = & (*sample_iter).second;
            grid_samples_vec.push_back(ws);
        }
        

        for (size_t s = 0; s != grid_samples_vec.size(); ++s)
        {
            fprintf(stdout, ":");
            WeightedSample & ws = * grid_samples_vec[s];
            for (size_t d = 0; d != ws.ndim; ++d)
            {
                fprintf(stdout, "\t%.10g\t%.10g", ws.x[d], ws.cdf[d]);
            }
            fprintf(stdout, "\n");
        }
    }


    return 0;
}
    
