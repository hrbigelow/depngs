struct output_dyn {
    struct output_dyn *next;
    char *out_buf; /* alloc'ed by thread */
    size_t out_size, out_alloc;
    int in_use;
};

struct thread_comp_input {
    posterior_wrapper *worker;
    char *in_buf;
    size_t in_size, in_alloc;

    // if any non-reference base has its test_quantile greater than
    // min_quantile_value it will be reported.
    float test_quantile, min_test_quantile_value;
    struct output_dyn *target;
};




int main()
{
    struct output_dyn *outdyn_pool, *outdyn_head;
    struct pthread_mutex_t *pool_mtx, *read_mtx;
    struct pthread_t *threads;
    struct thread_comp_input *worker_inputs;
    size_t num_threads, num_pool = num_threads * 2;

    worker_inputs = malloc(num_threads * sizeof(struct thread_comp_input));
    outdyn_pool = calloc(num_pool, sizeof(struct output_dyn));
    outdyn_head = NULL;
    threads = malloc(num_threads * sizeof(struct pthread_t));

    pthread_mutex_init(pool_mtx, mtx_attr);
    pthread_mutex_init(read_mtx, mtx_attr);

    size_t t, p;
    int rc;

    /* only need initialize the pointer to NULL */
    for (p = 0; p != num_threads + num_extra; ++p)
        outdyn_pool[p].out_buf = NULL;

    for (t = 0; t != num_threads; ++t)
    {
        worker_inputs[t].worker = 0;
        worker_inputs[t].test_quantile = 0;
        worker_inputs[t].min_test_quantile_value = 0;
        worker_inputs[t].target = NULL;

        rc = pthread_create(&threads[t], NULL, worker, (void *)&worker_inputs[t]);
        assert(rc == 0);
    }

    for (t = 0; t != num_threads; ++t)
    {
        rc = pthread_join(threads[t], NULL);
        assert(rc == 0);
    }

    
    pthread_mutex_destroy(pool_mtx);
    pthread_mutex_destroy(read_mtx);
    free(outdyn_pool);
    free(threads);
}

void *comp_worker(void *args)
{
    struct worker_input *par = args;
    double *sample_points_buf;

    while (1)
    {
        int rc;
        /* don't stop until you find a free dyn buffer */
        rc = pthread_mutex_lock(pool_mtx);
        while (! par->target)
            for (p = 0; p != num_pool; ++p)
                if (! outdyn_pool[p].in_use)
                {
                    par->target = &outdyn_pool[p];
                    break;
                }

        /* update the linked list */
        par->target.in_use = 1;
        par->target.next = outdyn_head;
        outdyn_head = par->target;
        rc = pthread_mutex_unlock(pool_mtx);

        /* read chunk into input buffer */
        rc = pthread_mutex_lock(read_mtx);
        read_chunk(ix, pileup_input_fh, 

        rc = pthread_mutex_unlock(read_mtx);

        if (1 /* no content read */)
            pthread_exit(NULL);

        sample_points_buf = 
            malloc(par->worker->s.final_num_points * 4 * sizeof(double));
        
        char 
            *line = par->in_buf, 
            *end = line + par->target->in_size, 
            *next;

        size_t max_output_line = 1000;
        while (line != end)
        {
            next = strchr(line, '\n') + 1;
            next[-1] = '\0'; /* weird, but works */

            ALLOC_GROW_REMAP(par->target->out_buf,
                             write_ptr,
                             write_ptr - par->target->out_buf + max_output_line,
                             par->target->out_alloc);

            write_ptr =
                par->worker->process_line_comp(line, write_ptr,
                                               sample_points_buf,
                                               par->test_quantile,
                                               par->min_test_quantile_value);
            line = next;
        
        }
        par->target->out_size = write_ptr - par->target->out_buf;
        par->target->in_use = 0;
        par->target = NULL;
    }
    free(sample_points_buf);
    
}
