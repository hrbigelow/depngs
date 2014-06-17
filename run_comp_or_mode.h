#ifndef _RUN_COMP_OR_MODE_H
#define _RUN_COMP_OR_MODE_H

#include <cstddef>

int run_comp_or_mode(size_t max_mem,
                     size_t num_threads,
                     size_t min_quality_score,
                     char const* label_string,
                     char const* quantiles_file,
                     char const* prior_alphas_file,
                     char const* pileup_input_file,
                     char const* jpd_data_params_file,
                     char const* posterior_output_file,
                     char const* cdfs_output_file,
                     size_t tuning_num_points,
                     size_t final_num_points,
                     bool verbose,
                     void * (*worker)(void *));

#endif // _RUN_COMP_OR_MODE_H
