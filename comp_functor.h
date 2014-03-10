#ifndef _COMP_FUNCTOR_H
#define _COMP_FUNCTOR_H

#include <vector>
#include <fstream>

class NucleotideStats;
class ErrorEstimate;
class Dirichlet;
class Posterior;
class Metropolis;


struct posterior_wrapper
{
    double mode_tolerance;
    size_t max_modefinding_iterations;
    size_t max_tuning_iterations;
    size_t tuning_num_points;
    size_t final_num_points;
    double initial_point[4];
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
    double * sample_points_buf;
    std::vector<double *> sample_points_sortable;

    char * label_string;

    NucleotideStats * params;
    ErrorEstimate * model;
    Dirichlet * prior;
    Posterior * posterior;
    Metropolis * sampler;

    posterior_wrapper(char const* jpd_data_params_file,
                      double * prior_alphas,
                      size_t min_quality_score,
                      double * quantiles,
                      size_t num_quantiles,
                      char const* label_string,
                      FILE * cdfs_output_fh,
                      pthread_mutex_t * file_writing_mutex,
                      bool compute_anomaly);

    ~posterior_wrapper();

    void process_line_comp(char const* pileup_line, char * out_buffer);
    void process_line_mode(char const* pileup_line, char * out_buffer);
};


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
