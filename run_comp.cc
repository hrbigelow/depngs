#include <cstddef>

#include "run_comp_or_mode.h"
#include "pileup_tools.h"
#include "comp_functor.h"
#include "file_utils.h"
#include "defs.h"
#include "cache.h"
#include "file_binary_search.h"

#include <gsl/gsl_errno.h>

int run_comp_or_mode(size_t max_mem,
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
    {
        quantiles = ParseNumbersFile(quantiles_file, & num_quantiles);
    }

    double prior_alphas[NUM_NUCS];
    for (size_t i = 0; i != NUM_NUCS; ++i)
        prior_alphas[i] = prior_alpha;
    
    

    FILE * posterior_output_fh = open_if_present(posterior_output_file, "w");
    FILE * cdfs_output_fh = open_if_present(cdfs_output_file, "w");
    FILE * pileup_input_fh = open_if_present(pileup_input_file, "r");

    size_t base_chunk_size = max_mem; /* user provided.  Must be >= 1e8 (100MB) */
    size_t max_pileup_line_size; /* defaulted to 10000.  automatically updated */

    // 0. parse contig_order file
    FILE *contig_order_fh = fopen(contig_order_file, "r");
    if (! contig_order_fh)
    {
        fprintf(stderr, "Couldn't open contig order file %s\n", contig_order_file);
        exit(1);
    }

    char contig[1024];
    long cix;
    unsigned index;
    while (! feof(contig_order_fh))
    {
        fscanf(contig_order_fh, "%s\t%u\n", contig, &index);
        dict_add_item(contig, index);
    }
    fclose(contig_order_fh);
    
    dict_build();

    /* 1. parse all query ranges into 'queries' and sort them */
    unsigned num_queries = 0, num_alloc = 10;
    struct locus_range *queries = 
        (struct locus_range *)malloc(num_alloc * sizeof(struct locus_range)),
        *qend,
        *q;
    
    /* construct the set of non-overlapping query ranges */
    char reformat_buf[1000];
    unsigned beg_pos, end_pos;
    while (fscanf(locus_fh, "%s\t%u\t%u\n", contig, &beg_pos, &end_pos) == 3)
    {
        sprintf(reformat_buf, "%s\t%u\t", contig, beg_pos);
        queries[num_queries].beg = init_locus(reformat_buf);

        sprintf(reformat_buf, "%s\t%u\t", contig, end_pos);
        queries[num_queries].end = init_locus(reformat_buf);

        ++num_queries;
        ALLOC_GROW(queries, num_queries + 1, num_alloc);
    }   
    fclose(locus_fh);

    qsort(queries, num_queries, sizeof(queries[0]), less_locus_range);
    
    /* 2. edit queries to eliminate interval overlap */
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
    size_t output_unit_size = 
        (strlen(label_string) + 112 + (10 * num_quantiles) + (14 + num_quantiles)) * 5;

    // the functions in the workers require this.
    gsl_set_error_handler_off();

    /* 3. create and initialize a root index node representing the
       entire pileup file */
    struct file_bsearch_index
        *root = find_root_index(pileup_input_fh),
        *ix = root;

    int new_query = 1;

    /* cast to ptrdiff_t so we can do signed-comparison, even though these
       are always positive. */
#define CHUNK_SIZE() (ptrdiff_t)(base_chunk_size + max_pileup_line_size)
#define FILE_SPAN() (ptrdiff_t)(end_off - start_off)
    
    /* already a ptrdiff_t type */
#define BASE_LEFT() (base_end - write_ptr)

    /* allocate a reasonable size for the output buffer */
    size_t t;
    size_t out_buf_size = 1e7;
    for (t = 0; t != num_threads; ++t)
    {
        worker_inputs[t].out_buf = (char *)malloc(out_buf_size);
        worker_inputs[t].out_alloc = out_buf_size;
        worker_inputs[t].out_size = 0;
    }

    /* redo the input strategy assuming off_index input */
    while (q != qend)
    {
        while (BASE_LEFT() > 0)
        {
            /* find file offsets for current query, creating index nodes and
               updating ix in the process */
            if (new_query)
            {
                ix = find_loose_index(ix, q->beg, pileup_input_fh);
                start_off = off_lower_bound(ix, q->beg);
                ix = find_loose_index(ix, q->end, pileup_input_fh);
                end_off = off_upper_bound(ix, q->end);

                fseeko(pileup_input_fh, start_off, SEEK_SET);
            }

            /* fill the buffer as much as possible with the next query range.
               afterwards, q now points to the next range to retrieve. */
            if (FILE_SPAN() < BASE_LEFT())
            {
                write_ptr += fread(write_ptr, 1, FILE_SPAN(), pileup_input_fh);
                ++q;
                new_query = 1;
            }
            else
            {
                /* read up to the base buffer, then read the next line
                   fragment.  realloc both chunk_buf and line_buf as
                   necessary */
                write_ptr += fread(write_ptr, 1, BASE_LEFT(), pileup_input_fh);
                size_t oldmax = max_pileup_line_size;
                ssize_t line_length = getline(&line_buf, &max_pileup_line_size, pileup_input_fh);
                assert(line_length >= 0);

                if (oldmax != max_pileup_line_size)
                {
                    /* should happen very rarely.  if so, realloc and
                       update buffers. */
                    size_t write_pos = write_ptr - chunk_buf;
                    chunk_buf = (char *)realloc(chunk_buf, CHUNK_SIZE());
                    base_end = chunk_buf + base_chunk_size;
                    write_ptr = chunk_buf + write_pos;
                }
                strcpy(write_ptr, line_buf);
                write_ptr += line_length;

                /* update q->beg to the position just after the last line read */
                char *last_line = (char *)memrchr(chunk_buf, '\n', write_ptr - chunk_buf - 1) + 1;
                q->beg = init_locus(prev);
                ++q->beg.pos;
                new_query = 0;
            }
        }

        /* at this point, the range [chunk_buf, write_ptr) contains a full
           set of lines, and q identifies the next set of loci to
           process.  process the loci in the range. */
        pthread_t *threads = (pthread_t *)malloc(sizeof(pthread_t) * num_threads);
        char *cut;
        
        worker_inputs[0].beg = read_beg;
        for (size_t t = 1; t != num_threads; ++t)
        {
            cut = chunk_buf + (write_ptr - chunk_buf) * t / num_threads;
            cut = strchr(cut, '\n') + 1;

            worker_inputs[t].beg = worker_inputs[t-1].end = cut;
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

        for (size_t t = 0; t < num_threads; ++t) {
            int rc = pthread_join(threads[t], NULL);
            assert(0 == rc);
        }

        free(threads);
        
        for (t = 0; t != num_threads; ++t)
            fwrite(worker_inputs[t].out_buf, 1, worker_inputs.out_size, posterior_output_fh);

        fflush(posterior_output_fh);


    }
    
    /* 
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
    */

    fclose(pileup_input_fh);
    free(chunk_buf);

    fclose(posterior_output_fh);
    if (cdfs_output_fh != NULL)
        fclose(cdfs_output_fh);

    delete quantiles;
    free(worker_inputs);

    return 0;
    
}
