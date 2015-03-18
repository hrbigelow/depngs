#ifndef _RUN_COMP_H
#define _RUN_COMP_H

#include <cstddef>

struct posterior_settings;

int run_comp(size_t max_mem,
             size_t num_threads,
             const char *fastq_type,
             const char *label_string,
             const char *quantiles_file,
             const char *pileup_input_file,
             const char *contig_order_file,
             const char *query_range_file,
             const char *jpd_data_params_file,
             const char *posterior_output_file,
             const char *cdfs_output_file,
             posterior_settings *pset,
             bool verbose);

#endif // _RUN_COMP_H
