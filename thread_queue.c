/* managed output buffers */
struct output_node {
    struct output_node *next;
    char *out_buf; /* alloc'ed by thread */
    size_t out_size, out_alloc;
    int in_use;
};


struct thread_comp_input {
    void *par; /* generic parameters to be passed to the thread worker */
    char *in_buf;
    size_t in_size, in_alloc;
    struct output_node *target;
    struct thread_queue *tq;
};


/* the main configuration object. */
struct thread_queue {
    void *(*worker)(void *); /* generic client-provided function for processing*/
    void (*reader)(FILE *); /* */
    struct output_node *out_pool, *out_head;
    struct pthread_mutex_t *pool_mtx, *read_mtx;
    struct pthread_t *threads;
    size_t num_threads, num_pool;
};


/* initialize resources */
struct thread_queue *
thread_queue_init(void *(*worker)(void *),
                  void (*reader)(FILE *fh),
                  size_t num_threads,
                  size_t num_extra_in_pool)
{
    worker_inputs = malloc(num_threads * sizeof(struct thread_comp_input));
    out_pool = calloc(num_pool, sizeof(struct output_node));
    out_head = NULL;
    threads = malloc(num_threads * sizeof(struct pthread_t));

    pthread_mutex_init(pool_mtx, mtx_attr);
    pthread_mutex_init(read_mtx, mtx_attr);

    size_t t, p;
    int rc;

    /* only need initialize the pointer to NULL */
    for (p = 0; p != num_threads + num_extra; ++p)
        out_pool[p].out_buf = NULL;

}


/* start things running. */
int thread_queue_run(struct thread_queue *tq)
{

    size_t t;
    int rc;
    for (t = 0; t != tq->num_threads; ++t) {
        rc = pthread_create(&tq->threads[t], NULL, 
                            worker_func,
                            (void *)&worker_inputs[t]);
        if (rc)
        {
            fprintf(stderr, "Error: pthread_create returned value %i at %s:%u\n", 
                    rc, __FILE__, __LINE__);
            exit(1);
        }

        rc = pthread_join(tq->threads[t], NULL);
        if (rc)
        {
            fprintf(stderr, "Error: pthread_join returned value %i at %s:%u\n", 
                    rc, __FILE__, __LINE__);
            exit(1);
        }
    }
    return 0; /* not sure what else to return */
}


/* release resources */
int thread_queue_free(struct thread_queue *tq)
{
    pthread_mutex_destroy(tq->pool_mtx);
    pthread_mutex_destroy(tq->read_mtx);
    free(tq->out_pool);
    free(tq->threads);
}


static void *worker_func(void *args)
{
    struct thread_comp_input *par = args;
    struct thread_queue *tq = par->tq;
    size_t p;
    while (1)
    {
        int rc;
        /* don't stop until you find a free dyn buffer */
        rc = pthread_mutex_lock(tq->pool_mtx);
        while (! par->target)
            for (p = 0; p != num_pool; ++p)
                if (! tq->out_pool[p].in_use)
                {
                    par->target = &tq->out_pool[p];
                    break;
                }

        /* update the linked list */
        par->target.in_use = 1;
        par->target.next = tq->out_head;
        tq->out_head = par->target;
        rc = pthread_mutex_unlock(tq->pool_mtx);

        /* read chunk into input buffer */
        rc = pthread_mutex_lock(tq->read_mtx);
        tq->reader();

        rc = pthread_mutex_unlock(tq->read_mtx);

        if (1 /* no content read */)
            pthread_exit(NULL);

        tq->worker();

        /* */

        par->target->out_size = write_ptr - par->target->out_buf;
        par->target->in_use = 0;
        par->target = NULL;
    }
    free(sample_points_buf);
    
}
