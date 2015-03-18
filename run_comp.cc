#include "run_comp.h"
#include "pileup_tools.h"
#include "comp_worker.h"
#include "yepLibrary.h"

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
             const char *fastq_type,
             const char *label_string,
             const char *quantiles_file,
             const char *pileup_input_file,
             const char *contig_order_file,
             const char *query_range_file,
             const char *jpd_data_params_file,
             const char *posterior_output_file,
             const char *cdfs_output_file,
             struct posterior_settings *pset,
             bool verbose)
{
    double *quantiles;
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

    FILE *posterior_output_fh = open_if_present(posterior_output_file, "w");
    FILE *cdfs_output_fh = open_if_present(cdfs_output_file, "w");
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

    struct sample_attributes sample_atts;
    init_sample_attributes(jpd_data_params_file,
                           label_string,
                           pileup_input_file,
                           &sample_atts);

    /* 1. create and initialize a root index node representing the
       entire pileup file */
    size_t scan_thresh_size = 1e6;
    file_bsearch_init(init_locus, scan_thresh_size);
    
    struct file_bsearch_index ix = file_bsearch_make_index(sample_atts.fh);

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
        queries[0].end.lo--; /* necessary for this pseudo-query to fit in the index root */
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

    // create the reusable resources here
    struct comp_worker_input *worker_buf =
        (struct comp_worker_input *)
        malloc(n_threads * sizeof(struct comp_worker_input));

    void **worker_inputs = (void **)malloc(n_threads * sizeof(void *));

    size_t t;
    for (t = 0; t != n_threads; ++t)
    {
        worker_buf[t].sample_atts = sample_atts;
        worker_buf[t].pset = *pset;
        memcpy(worker_buf[t].quantiles, quantiles, sizeof(double) * n_quantiles);
        worker_buf[t].n_quantiles = n_quantiles;
        // worker_buf[t].test_quantile = test_quantile;
        // worker_buf[t].min_test_quantile_value = min_test_quantile_value;
    }

    for (t = 0; t != n_threads; ++t)
        worker_inputs[t] = &worker_buf[t];
    
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

    enum YepStatus status = yepLibrary_Init();
    assert(status == YepStatusOk);

    thread_queue_run(tqueue);

    thread_queue_free(tqueue);

    fclose(sample_atts.fh);
    fclose(posterior_output_fh);
    if (cdfs_output_fh != NULL)
        fclose(cdfs_output_fh);

    delete quantiles;
    free(worker_buf);
    free(worker_inputs);

    return 0;
    
}
