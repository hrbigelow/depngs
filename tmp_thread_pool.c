
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

        par->target->out_size = write_ptr - par->target->out_buf;
        par->target->in_use = 0;
        par->target = NULL;
    }
    free(sample_points_buf);
    
}



    /* 
       Responsibilities of the client-provided worker:

       char *worker(void *param, const char *in_beg, const char *in_end,
                    char **out_buf, size_t *out_alloc);

       1.  Given out_buf, out_alloc
       2.  Write to out_buf, re-alloc and update out_alloc as necessary
       3.  Return pointer inside out_buf signifying end of write position
       4.  Be provided input in the form of a range [beg, end) (should this be read-only?)
       5.  Accept a void * for any required parameters to do this job.

    */

    /* 
       Responsibilities of the client-provided reader:

       reader(char **in_buf, 
       1.  Should we assume a FILE?  For the moment, yes

     */

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







the generic worker function will be library-provided, and it will
be passed to the thread.  it thus must conform to the prototype of

void *worker(void *args)



Each thread, when created, needs access to the main configuration
object, plus its own input buffer.  



what are the inputs and outputs of the worker.  conceptually, the
worker needs to be given just an input range and tells us how
much it has output.

The reader also must be in charge of re-allocating the buffer,
since the decision to reallocate must be based on the structure
of the input itself.



