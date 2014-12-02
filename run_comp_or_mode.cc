#include <cstddef>

#include "run_comp_or_mode.h"
#include "pileup_tools.h"
#include "comp_functor.h"
#include "file_utils.h"
#include "defs.h"

#include <gsl/gsl_errno.h>

int run_comp_or_mode(size_t max_mem,
                     size_t num_threads,
                     size_t min_quality_score,
                     char const* fastq_type,
                     char const* label_string,
                     char const* quantiles_file,
                     float prior_alpha,
                     char const* pileup_input_file,
                     char const* jpd_data_params_file,
                     char const* posterior_output_file,
                     char const* cdfs_output_file,
                     struct posterior_settings *pset,
                     double test_quantile,
                     double min_test_quantile_value,
                     bool verbose,
                     void * (*worker)(void *))
{
    double * quantiles;
    size_t num_quantiles;

    if (strcmp(quantiles_file, "/dev/null") == 0)
    {
        quantiles = new double[5];
        quantiles[0] = 0.005;
        quantiles[1] = 0.05;
        quantiles[2] = 0.5;
        quantiles[3] = 0.95;
        quantiles[4] = 0.995;
        num_quantiles = 5;
    }
    else
    {
        quantiles = ParseNumbersFile(quantiles_file, & num_quantiles);
    }

    double prior_alphas[NUM_NUCS];
    for (size_t i = 0; i != NUM_NUCS; ++i)
        prior_alphas[i] = prior_alpha;
    
    

    FILE * posterior_output_fh = open_if_present(posterior_output_file, "w");
    FILE * cdfs_output_fh = open_if_present(cdfs_output_file, "w");
    FILE * pileup_input_fh = open_if_present(pileup_input_file, "r");

    size_t chunk_size = max_mem;
    char * chunk_buffer_in = new char[chunk_size + 1];
    int offset;

    if (fastq_type)
        offset = fastq_type_to_offset(fastq_type);
    else
        offset = fastq_offset(pileup_input_file, chunk_buffer_in, chunk_size);

    if (offset == -1)
    {
        fprintf(stderr, "Could not determine fastq type of this pileup file.\n");
        return 1;
    }

    PileupSummary::set_offset(offset);

    //we are integrating the actual posterior
    // char const* dimension_labels[] = { "A", "C", "G", "T" };
    // char line_label[1000];

    // double mode_tolerance = 1e-60;
    // size_t max_modefinding_iterations = 3000;
    // double initial_point[] = { 0.25, 0.25, 0.25, 0.25 };

    // Create strand-marginal model parameters for use with the anomaly scoring
    // double * pos_strand_marginal_params = new double[Nucleotide::num_bqs * 4];
    // double * neg_strand_marginal_params = new double[Nucleotide::num_bqs * 4];

    size_t nbytes_read, nbytes_unused = 0;
    char * last_fragment;
    char * read_pointer = chunk_buffer_in;


    // this will be allocated / deallocated once for each chunk processed
    char * chunk_buffer_out;
    char ** output_lines;

    // size_t initial_autocor_offset = 30;
    
    // used for slice sampling
    // size_t num_bits_per_dim = 62;
    // bool is_log_integrand = true;
    // size_t initial_sampling_range = 62 * truncated_ndim;
    // bool may_underflow = true;
    // double prior_alpha0 = std::accumulate(prior_alphas, prior_alphas + 4, 0.0);
    // bool use_independence_chain_mh = true;

    // create the reusable resources here
    wrapper_input *worker_inputs = new wrapper_input[num_threads];

    pthread_mutex_t file_writing_mutex;

    for (size_t t = 0; t != num_threads; ++t)
    {
        worker_inputs[t].worker = 
            new posterior_wrapper(jpd_data_params_file,
                                  prior_alphas,
                                  min_quality_score,
                                  quantiles,
                                  num_quantiles,
                                  label_string,
                                  cdfs_output_fh,
                                  & file_writing_mutex,
                                  *pset,
                                  verbose);
    }
    
    // number of characters needed for a single inference
    size_t output_unit_size = 
        (strlen(label_string) + 112 + (10 * num_quantiles) + (14 + num_quantiles)) * 5;

    // the functions in the workers require this.
    gsl_set_error_handler_off();

    /* redo the input strategy assuming off_index input */
    

    while (! feof(pileup_input_fh))
    {
        nbytes_read = fread(read_pointer, 1, chunk_size - nbytes_unused, pileup_input_fh);

        std::vector<char *> pileup_lines =
            FileUtils::find_complete_lines_nullify(chunk_buffer_in, &last_fragment);

        chunk_buffer_out = new char[output_unit_size * pileup_lines.size() + 1];
        output_lines = new char*[pileup_lines.size()];

        char * current_line = chunk_buffer_out;
        for (size_t l = 0; l != pileup_lines.size(); ++l)
        {
            output_lines[l] = current_line;
            current_line += output_unit_size;
        }

        read_pointer[nbytes_read] = '\0';

        // here, create N pthreads.  Each should be
        pthread_t *threads = new pthread_t[num_threads];
        size_t worker_load = pileup_lines.size() / num_threads;
        for (size_t t = 0; t != num_threads; ++t)
        {
            worker_inputs[t].beg = pileup_lines.begin() + (t * worker_load);
            worker_inputs[t].end = (t == num_threads - 1) 
                ? pileup_lines.end()
                : pileup_lines.begin() + ((t + 1) * worker_load);
            worker_inputs[t].out_start = output_lines + (t * worker_load);
            worker_inputs[t].test_quantile = test_quantile;
            worker_inputs[t].min_test_quantile_value = min_test_quantile_value;

            int rc = pthread_create(&threads[t], NULL, 
                                    worker,
                                    static_cast<void *>(& worker_inputs[t]));
            assert(rc == 0);
        }

        for (size_t t = 0; t < num_threads; ++t) {
            int rc = pthread_join(threads[t], NULL);
            assert(0 == rc);
        }
        
        // write the buffers
        for (size_t l = 0; l != pileup_lines.size(); ++l)
            fwrite(output_lines[l], 1, strlen(output_lines[l]), posterior_output_fh);

        fflush(posterior_output_fh);

        nbytes_unused = strlen(last_fragment);
        memmove(chunk_buffer_in, last_fragment, nbytes_unused);
        read_pointer = chunk_buffer_in + nbytes_unused;
        delete chunk_buffer_out;
        delete output_lines;
        delete threads;

    }
    fclose(pileup_input_fh);
    delete chunk_buffer_in;

    fclose(posterior_output_fh);
    if (cdfs_output_fh != NULL)
    {
        fclose(cdfs_output_fh);
    }

    delete quantiles;
    delete worker_inputs;

    return 0;
    
}
