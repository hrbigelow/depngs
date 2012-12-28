#include <cmath>
#include <ctime>
#include <cstdio>

#include <gsl/gsl_statistics_double.h>

#include "henry/stats_tools.h"


double distance(double const* a, double const* b, size_t num_elem)
{
    double sum_sq = 0.0;
    for (size_t i = 0; i != num_elem; ++i)
    {
        sum_sq += (a[i] - b[i])*(a[i] - b[i]);
    }
    return sqrt(sum_sq);
}


int main()
{
    srand(time(NULL));
    size_t num_x = 100;
    double * x = new double[num_x];
    double * y = new double[num_x];

    double inc_mean = 0.0;
    double mean;
    double dif;

    printf("adding to means\n");
    for (size_t i = 0; i != num_x; ++i)
    {
        x[i] = static_cast<double>(rand()) / static_cast<double>(RAND_MAX);
        y[i] = static_cast<double>(rand()) / static_cast<double>(RAND_MAX);
        inc_mean = add_sample_to_mean(inc_mean, i, x[i]);
        mean = gsl_stats_mean(x, 1, i + 1);
        dif = inc_mean - mean;
        printf("mean diff: %Zu\t%g\n", i, dif);
    }

    
    printf("removing from means\n");
    for (size_t i = 0; i != num_x - 1; ++i)
    {
        size_t remain = num_x - (i + 1);
        inc_mean = remove_sample_from_mean(inc_mean, remain + 1, x[i]);
        mean = gsl_stats_mean(x + (i+1), 1, remain);
        dif = inc_mean - mean;
        printf("mean diff: %Zu\t%g\n", i, dif);
    }

        
    printf("adding to covariances\n");

    double covariance;
    double mean_x = gsl_stats_mean(x, 1, 2);
    double mean_y = gsl_stats_mean(y, 1, 2);
    double inc_covariance = gsl_stats_covariance_m(x, 1, y, 1, 2, mean_x, mean_y);
    
    for (size_t i = 2; i != num_x; ++i)
    {
        inc_covariance = 
            add_sample_to_covariance(inc_covariance, i,
                                     x[i], mean_x, y[i], mean_y);

        mean_x = add_sample_to_mean(mean_x, i, x[i]);
        mean_y = add_sample_to_mean(mean_y, i, y[i]);

        covariance = gsl_stats_covariance_m(x, 1, y, 1, i + 1, mean_x, mean_y);

        dif = inc_covariance - covariance;
        printf("covar diff: %Zu\t%g\n", i, dif);
    }    

    mean_x = gsl_stats_mean(x, 1, num_x);
    mean_y = gsl_stats_mean(y, 1, num_x);

    inc_covariance = 
        gsl_stats_covariance_m(x, 1, y, 1, num_x, mean_x, mean_y);

    printf("removing from covariances\n");
    for (size_t i = 0; i != num_x - 2; ++i)
    {
        size_t remain = num_x - (i + 1);

        inc_covariance = 
            remove_sample_from_covariance(inc_covariance, remain + 1,
                                          x[i], mean_x, y[i], mean_y);

        mean_x = remove_sample_from_mean(mean_x, remain + 1, x[i]);
        mean_y = remove_sample_from_mean(mean_y, remain + 1, y[i]);

        covariance = 
            gsl_stats_covariance_m(x + (i+1), 1, y + (i+1), 1, remain, mean_x, mean_y);

        dif = inc_covariance - covariance;
        printf("covar diff: %Zu\t%g\n", i, dif);
    }    

    printf("adding to means and covariances");
    double mean_vector_inc[2], mean_vector[2];
    double covariance_matrix_inc[4], covariance_matrix[4];
    double * x_flat = new double[num_x * 2];

    size_t ndim = 2;

    for (size_t i = 0; i != num_x * ndim; ++i)
    {
        x_flat[i] = static_cast<double>(rand()) / static_cast<double>(RAND_MAX);
    }   
 
    size_t num_initial_samples = 2;
    multivariate_mean_covariance(x_flat, ndim, num_initial_samples, 
                                 mean_vector_inc, covariance_matrix_inc);

    for (size_t i = 3; i != num_x; ++i)
    {
        multivariate_mean_covariance(x_flat, ndim, i, mean_vector, covariance_matrix);

        add_to_mean_covariance_matrix(x_flat + (ndim * (i-1)), ndim, i - 1, 
                                      mean_vector_inc, covariance_matrix_inc);
        double mean_dif = distance(mean_vector_inc, mean_vector, ndim);
        double covar_dif = distance(covariance_matrix_inc, covariance_matrix, ndim * ndim);
        printf("mean_covar diff: %Zu\t%g\t%g\n", i, mean_dif, covar_dif);
    }
    

    printf("removing from means and covariances");
    num_initial_samples = num_x;
    multivariate_mean_covariance(x_flat, ndim, num_initial_samples, 
                                 mean_vector_inc, covariance_matrix_inc);

    for (size_t i = 1; i != num_x - 3; ++i)
    {
        size_t remain = num_x - i;

        multivariate_mean_covariance(x_flat + (i * ndim), ndim, remain,
                                     mean_vector, covariance_matrix);
        
        remove_from_mean_covariance_matrix(x_flat + (ndim * (i - 1)), ndim, remain + 1, 
                                           mean_vector_inc, covariance_matrix_inc);

        double mean_dif = distance(mean_vector_inc, mean_vector, ndim);
        double covar_dif = distance(covariance_matrix_inc, covariance_matrix, ndim * ndim);
        printf("mean_covar diff: %Zu\t%g\t%g\n", i, mean_dif, covar_dif);
    }
    

    delete x;
    delete y;
    delete x_flat;
}
