#include <vector>
#include <utility>


#include "henry/tools.h"
#include "henry/error_estimate.h"
#include "henry/integrands.h"
#include "henry/metropolis.h"
#include "henry/sampling.h"
#include "henry/pileup_tools.h"
#include "henry/stats_tools.h"


int main(int argc, char ** argv)
{

    size_t tuning_num_points = static_cast<size_t>(atof(argv[1]));
    size_t final_num_points = static_cast<size_t>(atof(argv[2]));
    char * low_bounds_file = argv[3];
    char * hi_bounds_file = argv[4];
    char * quantiles_file = argv[5];

    printf ("tuning_num_points: %Zu\n"
            "final_num_points: %Zu\n"
            "low_bounds_file: %s\n"
            "hi_bounds_file: %s\n"
            "quantiles_file: %s\n",

            tuning_num_points,
            final_num_points,
            low_bounds_file,
            hi_bounds_file,
            quantiles_file
            );
    
    //parse mass fractions file
    double quantiles[1000];
    size_t num_quantiles = ParseNumbersFile(quantiles_file, quantiles);

    double low_bounds[1000];
    size_t num_dimensions = 
        ParseNumbersFile(low_bounds_file, low_bounds);

    double hi_bounds[1000];
    ParseNumbersFile(hi_bounds_file, hi_bounds);

    StepFunction step_function(low_bounds, hi_bounds, num_dimensions);

    double * norm_mean = new double[num_dimensions];
    std::fill(norm_mean, norm_mean + num_dimensions, 
              1.0 / static_cast<double>(num_dimensions));

    double * norm_covariance = new double[num_dimensions * num_dimensions];
    for (size_t d1 = 0; d1 != num_dimensions; ++d1)
    {
        for (size_t d2 = 0; d2 != num_dimensions; ++d2)
        {
            norm_covariance[d1 * num_dimensions + d2] = (d1 == d2) ? 1.0 : 0.0;
        }
    }

    bool is_log2_integrand = false;

    ProposalGaussianMH proposal(&step_function, norm_mean,
                                norm_covariance, num_dimensions,
                                is_log2_integrand);
        
    proposal.Init();

    size_t max_tuning_iterations = 100;
    size_t averaging_window_size = 10;
    size_t front_back_window_distance = 20;
    size_t autocor_max_offset = 10;
    size_t autocor_max_val = 0.1;

    Metropolis tuning_metropolis(num_dimensions, tuning_num_points);

    double * starting_point = new double[num_dimensions];
    for (size_t d = 0; d != num_dimensions; ++d)
    {
        starting_point[d] = (low_bounds[d] + hi_bounds[d]) / 2.0;
        //starting_point[d] = low_bounds[d];
    }

    tuning_metropolis.set_current_point(starting_point);

    tuning_metropolis.tune_proposal(&step_function, &proposal, 
                                    tuning_num_points,
                                    max_tuning_iterations,
                                    averaging_window_size,
                                    front_back_window_distance,
                                    autocor_max_offset);

    size_t best_autocor_offset =
        best_autocorrelation_offset(tuning_metropolis.get_samples(),
                                    num_dimensions, tuning_num_points,
                                    autocor_max_offset, autocor_max_val);
        
    size_t final_every_nth = best_autocor_offset;

    Metropolis final_metropolis(num_dimensions, final_num_points);

    final_metropolis.set_current_point(tuning_metropolis.get_current_point());

    double proposal_mean, proposal_variance;

    final_metropolis.sample(&step_function, &proposal,
                            final_num_points, 0, final_every_nth,
                            &proposal_mean, &proposal_variance);
        

    double * current_mean = new double[num_dimensions];
    double * current_covariance = new double[num_dimensions * num_dimensions];

    multivariate_mean_covariance(final_metropolis.get_samples(),
                                 num_dimensions, final_num_points,
                                 current_mean, current_covariance);

    fprintf(stdout, "Final statistics of sample points\n");
    print_mean_covariance(stdout, current_mean, current_covariance, 
                          num_dimensions);

    fprintf(stdout, "\n\n");
    fflush(stdout);

    std::vector<double *> final_points_sortable(final_num_points);
    for (size_t i = 0; i != final_num_points; ++i)
    {
        final_points_sortable[i] = 
            final_metropolis.get_samples() + (num_dimensions * i);
    }

    print_cdf_comparison(stdout, &step_function,
                         &final_points_sortable,
                         quantiles, num_quantiles, 
                         num_dimensions);

    return 0;

}
