#ifndef _COMP_FUNCTOR_H
#define _COMP_FUNCTOR_H

#include <vector>
#include <fstream>

class NucleotideStats;
class ErrorEstimate;
class Dirichlet;
class Posterior;
class Metropolis;
class SliceSampling;
class PileupSummary;

struct posterior_wrapper
{
    double mode_tolerance;
    size_t max_modefinding_iterations;
    size_t max_tuning_iterations;
    size_t tuning_num_points;
    size_t final_num_points;
    double initial_point[4];
    double mode_point[4];
    bool on_zero_boundary[4];
    double autocor_max_offset;
    double autocor_max_value;
    size_t initial_autocor_offset;
    size_t target_autocor_offset;
    size_t num_bits_per_dim;
    bool is_log_integrand;
    size_t initial_sampling_range;
    bool may_underflow;
    bool use_independence_chain_mh;
    bool verbose;
    double prior_alpha0;
    size_t min_quality_score;
    size_t num_quantiles;
    double * quantiles;
    FILE * cdfs_output_fh;
    pthread_mutex_t * file_writing_mutex;
    bool compute_anomaly;
    std::vector<double *> sample_points_sortable;

    char * label_string;

    NucleotideStats * params;
    ErrorEstimate * model;
    Dirichlet * prior;
    /* Posterior * posterior; */
    Metropolis * sampler;
    SliceSampling * slice_sampler;

    posterior_wrapper(char const* jpd_data_params_file,
                      double * prior_alphas,
                      size_t min_quality_score,
                      double * quantiles,
                      size_t num_quantiles,
                      char const* label_string,
                      FILE * cdfs_output_fh,
                      pthread_mutex_t * file_writing_mutex,
                      size_t tuning_num_points,
                      size_t final_num_points,
                      bool verbose);

    ~posterior_wrapper();

    size_t tune_mh(PileupSummary * locus, double * sample_points_buf);
    size_t tune_ss(double * sample_points_buf);
    void sample(PileupSummary * locus, double * sample_points_buf, char * algorithm_used);
    void values(double * points, size_t num_points, double * values);

    char * print_quantiles(PileupSummary * locus, 
                           char * algorithm_used, 
                           double * sample_points_buf,
                           char * out_buffer);


    char * process_line_comp(char const* pileup_line, char * out_buffer, double * sample_points_buf);
    char * process_line_mode(char const* pileup_line, char * out_buffer);
};


size_t discrete_comp_locus_bytes(size_t num_discrete_values);

char * print_discrete_comp(PileupSummary * locus,
                           char const* sample_label,
                           double * discrete_values,
                           size_t num_discrete_values,
                           size_t num_discrete_points_to_print,
                           double min_value_to_print,
                           char * out_buf);


struct wrapper_input
{
    posterior_wrapper * worker;
    std::vector<char *>::iterator beg;
    std::vector<char *>::iterator end;
    char ** out_start;
};


void * comp_worker(void * args);
void * mode_worker(void * args);

#endif // _COMP_FUNCTOR_H