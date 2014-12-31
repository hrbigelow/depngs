#include "run_comp.h"
#include "pileup_tools.h"
#include "comp_worker.h"

#include "defs.h"
#include "tools.h"

#include <string.h>

extern "C" {
#include "cache.h"
#include "dict.h"
#include "locus.h"
#include "range_line_reader.h"
}

#include <gsl/gsl_errno.h>


void result_offload(void *par, const struct managed_buf *buf)
{
    FILE **fh = (FILE **)par;
    fwrite(buf->buf, 1, buf->size, *fh);
}


int run_comp(size_t max_mem,
             size_t n_threads,
             size_t min_quality_score,
             const char *fastq_type,
             const char *label_string,
             const char *quantiles_file,
             float prior_alpha,
             const char *pileup_input_file,
             const char *contig_order_file,
             const char *query_range_file,
             const char *jpd_data_params_file,
             const char *posterior_output_file,
             const char *cdfs_output_file,
             struct posterior_settings *pset,
             double test_quantile,
             double min_test_quantile_value,
             bool verbose)
{
    double * quantiles;
    size_t n_quantiles;

    if (strcmp(quantiles_file, "/dev/null") == 0)
    {
        quantiles = new double[5];
        quantiles[0] = 0.005;
        quantiles[1] = 0.05;
        quantiles[2] = 0.5;
        quantiles[3] = 0.95;
        quantiles[4] = 0.995;
        n_quantiles = 5;
    }
    else
        quantiles = ParseNumbersFile(quantiles_file, & n_quantiles);

    double prior_alphas[NUM_NUCS];
    for (size_t i = 0; i != NUM_NUCS; ++i)
        prior_alphas[i] = prior_alpha;
    
    

    FILE *posterior_output_fh = open_if_present(posterior_output_file, "w");
    FILE *cdfs_output_fh = open_if_present(cdfs_output_file, "w");
    FILE *pileup_input_fh = open_if_present(pileup_input_file, "r");
    FILE *contig_order_fh = open_if_present(contig_order_file, "r");

    size_t base_chunk_size = max_mem; /* user provided.  Must be >= 1e8 (100MB) */
    
    /* 0. parse contig order file */
    char contig[1000];
    unsigned contig_index;
    while (! feof(contig_order_fh))
    {
        int n = fscanf(contig_order_fh, "%s\t%u\n", contig, &contig_index);
        if (n != 2)
        {
            fprintf(stderr, "Error: contig order file %s doesn't have the proper format\n",
                    contig_order_file);
            exit(1);
        }
        dict_add_item(contig, contig_index);
    }
    fclose(contig_order_fh);
    
    dict_build();


    /* 1. create and initialize a root index node representing the
       entire pileup file */
    size_t scan_thresh_size = 1e6;
    file_bsearch_init(init_locus, scan_thresh_size);
    
    struct file_bsearch_index ix = file_bsearch_make_index(pileup_input_fh);

    struct pair_ordering_range *queries, *q, *qend;
    size_t n_queries;

    if (query_range_file)
        queries = parse_query_ranges(query_range_file, &n_queries);
    else
    {
        /* simply set the 'query' to the entire file */
        n_queries = 1;
        queries = (struct pair_ordering_range *)
            malloc(sizeof(struct pair_ordering_range));
        queries[0] = { ix.root->span_beg, ix.root->span_end };
    }
    q = queries;
    qend = queries + n_queries;
    
    int offset;

    if (fastq_type)
        offset = fastq_type_to_offset(fastq_type);
    else
    {
        char *chunk_buf = (char *)malloc(base_chunk_size);
        offset = fastq_offset(pileup_input_file, chunk_buf, base_chunk_size);
        free(chunk_buf);
    }

    if (offset == -1)
    {
        fprintf(stderr, "Could not determine fastq type of this pileup file.\n");
        return 1;
    }

    PileupSummary::set_offset(offset);

    // size_t initial_autocor_offset = 30;
    
    // used for slice sampling
    // size_t num_bits_per_dim = 62;
    // bool is_log_integrand = true;
    // size_t initial_sampling_range = 62 * truncated_ndim;
    // bool may_underflow = true;
    // double prior_alpha0 = std::accumulate(prior_alphas, prior_alphas + 4, 0.0);
    // bool use_independence_chain_mh = true;

    // create the reusable resources here
    struct comp_worker_input **worker_inputs =
        (struct comp_worker_input **)
        malloc(n_threads * sizeof(struct comp_worker_input *));

    pthread_mutex_t file_writing_mutex;

    size_t t;
    for (t = 0; t != n_threads; ++t)
    {
        worker_inputs[t] = (struct comp_worker_input *)malloc(sizeof(struct comp_worker_input));
        worker_inputs[t]->worker = 
            new posterior_wrapper(jpd_data_params_file,
                                  prior_alphas,
                                  min_quality_score,
                                  quantiles,
                                  n_quantiles,
                                  label_string,
                                  cdfs_output_fh,
                                  & file_writing_mutex,
                                  *pset,
                                  verbose);
        worker_inputs[t]->sample_points_buf = 
            (double *)malloc(pset->final_n_points * 4 * sizeof(double));
        worker_inputs[t]->test_quantile = test_quantile;
        worker_inputs[t]->min_test_quantile_value = min_test_quantile_value;
    }
    
    // the functions in the workers require this.
    gsl_set_error_handler_off();
    
    struct range_line_reader_par reader_par = {
        &ix, 1, q, qend, init_locus, 1
    };

    /* theoretically, this will allow for one chunk to take twice as
       long as all the others, and still not cause a thread to wait
       for a free buffer. */
    size_t n_extra = n_threads;

    struct thread_queue *tqueue = 
        thread_queue_init(range_line_reader,
                          &reader_par,
                          comp_worker,
                          worker_inputs,
                          result_offload,
                          &posterior_output_fh,
                          n_threads,
                          n_extra,
                          1, /* number of input files */
                          1, /* number of output files */
                          max_mem);

    thread_queue_run(tqueue);

    thread_queue_free(tqueue);

    fclose(pileup_input_fh);
    fclose(posterior_output_fh);
    if (cdfs_output_fh != NULL)
        fclose(cdfs_output_fh);

    delete quantiles;
    for (t = 0; t != n_threads; ++t)
    {
        free(worker_inputs[t]->sample_points_buf);
        free(worker_inputs[t]);
    }

    free(worker_inputs);

    return 0;
    
}
