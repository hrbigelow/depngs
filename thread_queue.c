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


/* the main configuration object. */
struct thread_queue {
    void *(*worker)(void *); /* */
    void (*reader)(FILE *); /* */
    struct output_dyn *outdyn_pool, *outdyn_head;
    struct pthread_mutex_t *pool_mtx, *read_mtx;
    struct pthread_t *threads;
    size_t num_threads, num_pool;
};


/* initialize resources */
struct thread_queue *
thread_queue_init(void *(*worker)(void *),
                  void (*reader)(FILE *fh))
{
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

}


/* start things running. */
int thread_queue_run(struct thread_queue *tq)
{
}


/* release resources */
int thread_queue_free(struct thread_queue *tq)
{
    pthread_mutex_destroy(tq->pool_mtx);
    pthread_mutex_destroy(tq->read_mtx);
    free(tq->outdyn_pool);
    free(tq->threads);
}

