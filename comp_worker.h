#ifndef _COMP_FUNCTOR_H
#define _COMP_FUNCTOR_H

#include <vector>
#include <fstream>

extern "C" {
#include "thread_queue.h"
#include "ordering.h"
}

class NucleotideStats;
class ErrorEstimate;
class Dirichlet;
class Posterior;
class Metropolis;
class SliceSampling;
class PileupSummary;

enum sampling_method
    {
        METROPOLIS_HASTINGS = 'M',
        SLICE_SAMPLING = 'S',
        FAILED = '-'
    };

/* */
struct sample_details
{
    PileupSummary *locus;
    bool is_next;
    unsigned char dist_printed;
    double *sample_points;
    unsigned n_sample_points;
    sampling_method samp_method;
    bool mode_computed;
    size_t autocor_offset;
    char *current, *end;
    pair_ordering locus_ord;
    sample_details(void);
};


struct posterior_settings
{
    double gradient_tolerance;
    size_t max_modefinding_iterations;
    size_t max_tuning_iterations;
    size_t tuning_n_points;
    size_t final_n_points;
    double autocor_max_offset;
    double autocor_max_value;

    size_t initial_autocor_offset;
    size_t target_autocor_offset;
    bool is_log_integrand; // defaults to true
    size_t initial_sampling_range; // defaults to 62 * 3
};

// the mode that we use when there is no data
const double NULL_MODE[] = { 0.25, 0.25, 0.25, 0.25 };

struct posterior_wrapper
{
    struct posterior_settings s;
    double initial_point[4];
    double mode_point[4];
    bool on_zero_boundary[4];
    bool verbose;
    /* double prior_alpha0; */
    size_t min_quality_score;
    size_t n_quantiles;
    double *quantiles;
    FILE *cdfs_output_fh;
    pthread_mutex_t *file_writing_mutex;
    bool compute_anomaly;
    std::vector<double *> sample_points_sortable;

    char *label_string;

    NucleotideStats *params;
    ErrorEstimate *model;
    /* Dirichlet *prior; */
    Metropolis *sampler;
    SliceSampling *slice_sampler;

    posterior_wrapper(const char *jpd_data_params_file,
                      double *prior_alphas,
                      size_t min_quality_score,
                      double *quantiles,
                      size_t n_quantiles,
                      const char *label_string,
                      FILE *cdfs_output_fh,
                      pthread_mutex_t *file_writing_mutex,
                      posterior_settings s,
                      bool verbose);

    ~posterior_wrapper();

    void find_mode(void);
    void tune(sample_details *sd, double *estimated_mean);
    size_t tune_mh(PileupSummary *locus, double *sample_points_buf, double *estimated_mean);
    size_t tune_ss(PileupSummary *locus, double *sample_points_buf, double *estimated_mean);
    void sample(sample_details *sd, double *initial_point, size_t n_points);
    // void sample(PileupSummary *locus, double *sample_points_buf, char *algorithm_used);
    void values(double *points, size_t n_points, double *values);

    char *print_quantiles(sample_details *sd, char *out_buffer);

    char *process_line_comp(const char *pileup_line, char *out_buffer, 
                            double *sample_points_buf,
                            float test_quantile,
                            float min_test_quantile_value);
    
    char *process_line_mode(const char *pileup_line, char *out_buffer);
};


size_t discrete_comp_locus_bytes(size_t n_discrete_values);

char *print_discrete_comp(PileupSummary *locus,
                          const char *sample_label,
                          double *discrete_values,
                          size_t n_discrete_values,
                          size_t n_discrete_points_to_print,
                          double min_value_to_print,
                          char *out_buf);


struct comp_worker_input
{
    posterior_wrapper *worker;
    double *sample_points_buf; /* must be posterior_wrapper::final_n_points * 4 */

    /* if any non-reference base has its test_quantile greater than
       min_quantile_value it will be reported. */
    float test_quantile, min_test_quantile_value; 
};


thread_queue_worker_t comp_worker;

#endif // _COMP_FUNCTOR_H
