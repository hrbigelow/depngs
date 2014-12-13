#include <cstddef>

#include "run_comp.h"
#include "pileup_tools.h"
#include "comp_functor.h"
#include "file_utils.h"
#include "defs.h"
extern "C" {
#include "cache.h"
#include "dict.h"
#include "pileup_bsearch.h"
}

#include <gsl/gsl_errno.h>



int run_comp(size_t max_mem,
             size_t num_threads,
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
        quantiles = ParseNumbersFile(quantiles_file, & num_quantiles);

    double prior_alphas[NUM_NUCS];
    for (size_t i = 0; i != NUM_NUCS; ++i)
        prior_alphas[i] = prior_alpha;
    
    

    FILE *posterior_output_fh = open_if_present(posterior_output_file, "w");
    FILE *cdfs_output_fh = open_if_present(cdfs_output_file, "w");
    FILE *pileup_input_fh = open_if_present(pileup_input_file, "r");


    FILE *contig_order_fh = open_if_present(contig_order_file, "r");

    size_t base_chunk_size = max_mem; /* user provided.  Must be >= 1e8 (100MB) */
    size_t max_pileup_line_size; /* defaulted to 10000.  automatically updated */

    
    /* 0. parse contig order file */
    char contig[1024];
    unsigned index;
    while (! feof(contig_order_fh))
    {
        int n = fscanf(contig_order_fh, "%s\t%u\n", contig, &index);
        if (n != 2)
        {
            fprintf(stderr, "Error: contig order file %s doesn't have the proper format\n",
                    contig_order_file);
            exit(1);
        }
        dict_add_item(contig, index);
    }
    fclose(contig_order_fh);
    
    dict_build();

    struct locus_range {
        struct file_bsearch_ord beg, end;
    };

    struct locus_range *queries, *q, *qend;
    unsigned num_queries = 0, num_alloc;

    /* 1. create and initialize a root index node representing the
       entire pileup file */
    size_t scan_thresh_size = 1e6;
    file_bsearch_init(init_locus, scan_thresh_size);
    
    struct file_bsearch_index
        *root = find_root_index(pileup_input_fh),
        *ix = root;

    /* 2. parse all query ranges into 'queries' and sort them */
    char beg_contig[500], end_contig[500];
    if (query_range_file)
    {
        FILE *query_range_fh = open_if_present(query_range_file, "r");
        
        num_alloc = 10;
        queries = 
            (struct locus_range *)malloc(num_alloc * sizeof(struct locus_range));
        
        /* construct the set of non-overlapping query ranges */
        char reformat_buf[1000];
        unsigned beg_pos, end_pos;
        while (fscanf(query_range_fh, "%s\t%u\t%s\t%u\n", 
                      beg_contig, &beg_pos, end_contig, &end_pos) == 4)
        {
            sprintf(reformat_buf, "%s\t%u\t", beg_contig, beg_pos);
            queries[num_queries].beg = init_locus(reformat_buf);
            
            sprintf(reformat_buf, "%s\t%u\t", end_contig, end_pos);
            queries[num_queries].end = init_locus(reformat_buf);
            
            ++num_queries;
            ALLOC_GROW_TYPED(queries, num_queries + 1, num_alloc);
        }   
        fclose(query_range_fh);
    }        
    else
    {
        /* simply set the 'query' to the entire file */
        num_alloc = 1;
        num_queries = 1;
        queries = 
            (struct locus_range *)malloc(num_alloc * sizeof(struct locus_range));
        queries[0].beg = root->span.beg;
        queries[0].end = root->span.end;
        queries[0].end.lo--; /* so that root contains this query */
    }

    qsort(queries, num_queries, sizeof(queries[0]), less_locus_range);
    
    /* 3. edit queries to eliminate interval overlap */
    struct locus_range *p = NULL;
    for (q = queries; q != queries + num_queries - 1; ++q)
    {
        if (p && less_file_bsearch_ord(&p->end, &q->beg) > 0)
        {
            /* must be on same contig if they are overlapping and
               sorted */
            assert(p->end.hi == q->beg.hi);
            q->beg.lo = p->end.lo;
            if (q->end.lo < q->beg.lo)
                q->end.lo = q->beg.lo;
        }
        p = q;
    }

    q = queries;
    qend = queries + num_queries;
    
    char *chunk_buf = (char *)malloc(base_chunk_size),
        *base_end = chunk_buf + base_chunk_size,
        *write_ptr = chunk_buf,
        *line_buf = (char *)malloc(max_pileup_line_size);

    int offset;

    if (fastq_type)
        offset = fastq_type_to_offset(fastq_type);
    else
        offset = fastq_offset(pileup_input_file, chunk_buf, base_chunk_size);

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
    wrapper_input *worker_inputs = (wrapper_input *)malloc(num_threads * sizeof(wrapper_input));

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
    // size_t output_unit_size = 
    //     (strlen(label_string) + 112 + (10 * num_quantiles) + (14 + num_quantiles)) * 5;

    // the functions in the workers require this.
    gsl_set_error_handler_off();

    /* LEFT OFF */

    
        /* at this point, the range [chunk_buf, write_ptr) contains a full
           set of lines, and q identifies the next set of loci to
           process.  process the loci in the range. */
        pthread_t *threads = (pthread_t *)malloc(sizeof(pthread_t) * num_threads);
        char *cut;
        
        /* initialize ranges (notice the t = 1 initialization) */
        worker_inputs[0].beg = chunk_buf;
        for (size_t t = 1; t != num_threads; ++t)
        {
            cut = chunk_buf + (write_ptr - chunk_buf) * t / num_threads;
            cut = strchr(cut, '\n') + 1;
            worker_inputs[t].beg = worker_inputs[t-1].end = cut;
        }
        worker_inputs[num_threads - 1].end = write_ptr;

        /* initialize all other fieds (this time, t = 0 initialization) */
        for (size_t t = 0; t != num_threads; ++t)
        {
            worker_inputs[t].out_size = 0; /* we are re-freshing the output */

            /* out_buf and out_alloc are set before this, and are
               automatincally updated as needed during the worker
               execution. */
            worker_inputs[t].test_quantile = test_quantile;
            worker_inputs[t].min_test_quantile_value = min_test_quantile_value;

            int rc = pthread_create(&threads[t], NULL, 
                                    worker,
                                    static_cast<void *>(& worker_inputs[t]));
            assert(rc == 0);
        }

        for (size_t t = 0; t != num_threads; ++t) {
            int rc = pthread_join(threads[t], NULL);
            assert(rc == 0);
        }

        free(threads);
        
        for (t = 0; t != num_threads; ++t)
            fwrite(worker_inputs[t].out_buf, 1, worker_inputs[t].out_size, posterior_output_fh);

        fflush(posterior_output_fh);

        /* tells main loop we need to read more */
        write_ptr = chunk_buf;
    }
    
    fclose(pileup_input_fh);
    free(chunk_buf);

    fclose(posterior_output_fh);
    if (cdfs_output_fh != NULL)
        fclose(cdfs_output_fh);

    delete quantiles;
    free(worker_inputs);

    return 0;
    
}
