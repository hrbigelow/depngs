#include "thread_queue.h"

#include <pthread.h>
#include <stdlib.h>
#include <stdio.h>
#include <time.h>

enum buf_status {
    EMPTY,
    LOADING,
    FULL,
    UNLOADING
};

/* managed output buffers.  There will be N + E of these and they will
   be owned by 0 or 1 thread at any given time (but not necessarily
   the same thread each time. */
struct output_node {
    struct output_node *next;
    char *buf; /* alloc'ed by thread */
    size_t size, alloc;
    enum buf_status status;
};


/* input buffers.  There will be N of these, one for each thread, and
   each thread will use the same one for its lifespan. */
struct thread_comp_input {
    void *worker_par; /* generic parameters to be passed to the worker */
    char *buf;
    size_t size, alloc;
    struct output_node *out;
    struct thread_queue *tq;
};


/* the main configuration object. there will be just one of these for
   a given program (typically) or for a given task. */
struct thread_queue {
    thread_queue_reader_t *reader;
    thread_queue_worker_t *worker;
    thread_queue_offload_t *offload;

    /* parameters to be used by the reader within the read_mtx lock */
    void *reader_par, *offload_par;
    struct output_node *out_pool, *out_head, *out_tail;
    struct thread_comp_input *input;
    pthread_mutex_t pool_mtx, read_mtx;
    pthread_t *threads;
    /* struct pthread_attr_t thread_attr; */
    size_t num_threads, num_extra;
};


/* initialize resources */
struct thread_queue *
thread_queue_init(thread_queue_reader_t reader,
                  void *reader_par,
                  thread_queue_worker_t worker,
                  void *worker_par,
                  thread_queue_offload_t offload,
                  void *offload_par,
                  size_t num_threads,
                  size_t num_extra,
                  size_t max_input_mem)
{
    struct thread_queue *tq = malloc(sizeof(struct thread_queue));
    tq->reader = reader;
    tq->reader_par = reader_par;
    tq->worker = worker;
    tq->offload = offload;
    tq->offload_par = offload_par;
    tq->out_pool = calloc(num_threads + num_extra, sizeof(struct output_node));
    tq->out_head = NULL;
    tq->out_tail = NULL;
    tq->input = calloc(num_threads, sizeof(struct thread_comp_input));

    pthread_mutex_init(&tq->pool_mtx, NULL);
    pthread_mutex_init(&tq->read_mtx, NULL);
    tq->threads = malloc(num_threads * sizeof(pthread_t));
    tq->num_threads = num_threads;
    tq->num_extra = num_extra;


    /* only need initialize the pointer to NULL */
    size_t p, t;
    void **worker_pars = worker_par;
    for (t = 0; t != num_threads; ++t)
    {
        tq->input[t].worker_par = worker_pars[t];
        tq->input[t].alloc = max_input_mem / num_threads;
        tq->input[t].buf = malloc(tq->input[t].alloc);
        tq->input[t].size = 0;
        tq->input[t].tq = tq;
    }

    /* by default, use 2% of the memory of the input for the whole
       initial output */
    size_t out_chunk_size = max_input_mem * 2 / 100 / (num_threads + num_extra) + 10;
    for (p = 0; p != num_threads + num_extra; ++p)
    {
        tq->out_pool[p].alloc = out_chunk_size;
        tq->out_pool[p].size = 0;
        tq->out_pool[p].buf = malloc(tq->out_pool[p].alloc);
    }

    return tq;
}

static void *worker_func(void *args);


#define CHECK_THREAD(rc)                                                \
    if (rc) {                                                           \
        fprintf(stderr,                                                 \
                "Thread-related error at %s:%u with return code %i\n",  \
                __FILE__, __LINE__, rc);                                \
        exit(1);                                                        \
    }


/* start things running. */
int thread_queue_run(struct thread_queue *tq)
{

    size_t t;
    int rc;
    for (t = 0; t != tq->num_threads; ++t) {
        rc = pthread_create(&tq->threads[t], NULL, worker_func, &tq->input[t]);
        CHECK_THREAD(rc);
    }

    for (t = 0; t != tq->num_threads; ++t) {
        rc = pthread_join(tq->threads[t], NULL);
        CHECK_THREAD(rc);
    }
    return 0;
}


/* release resources */
void thread_queue_free(struct thread_queue *tq)
{
    pthread_mutex_destroy(&tq->pool_mtx);
    pthread_mutex_destroy(&tq->read_mtx);
    free(tq->threads);

    size_t t;
    for (t = 0; t != tq->num_threads; ++t)
        free(tq->input[t].buf);
    free(tq->input);

    for (t = 0; t != tq->num_threads + tq->num_extra; ++t)
        free(tq->out_pool[t].buf);
    free(tq->out_pool);
}


/* safely change the status flag of the 'out' node. */
void mutex_change_status(struct thread_queue *tq, struct output_node *out, enum buf_status status)
{
    int rc = pthread_mutex_lock(&tq->pool_mtx);
    CHECK_THREAD(rc);
    out->status = status;
    if (status == EMPTY)
        out->size = 0;

    rc = pthread_mutex_unlock(&tq->pool_mtx);
    CHECK_THREAD(rc);
}


#define ELAPSED_MS \
    ((((end_time).tv_sec * 1000000000 + (end_time).tv_nsec) -           \
      ((beg_time).tv_sec * 1000000000 + (beg_time).tv_nsec)) / 1000000)


#define PROGRESS_MSG(msg)                                   \
    do {                                                    \
        clock_gettime(CLOCK_REALTIME, &end_time);           \
        fprintf(stderr, "%li\t%li\t%s\n",                   \
                in - tq->input,                             \
                ELAPSED_MS,                                 \
                (msg));                                     \
        fflush(stderr);                                     \
    } while (0)

#define PROGRESS_START() clock_gettime(CLOCK_REALTIME, &beg_time)
    

/* this function is run by the thread */
static void *worker_func(void *args)
{
    struct thread_comp_input *in = args;
    struct thread_queue *tq = in->tq;
    size_t p, num_pool = tq->num_threads + tq->num_extra;

    struct timespec beg_time, end_time;

    while (1)
    {
        int rc;
        /* read chunk into input buffer.  this step can be done
           independent of any locking of the out-buffer pool. */
        rc = pthread_mutex_lock(&tq->read_mtx);
        CHECK_THREAD(rc);

        PROGRESS_START();
        tq->reader(tq->reader_par, &in->buf, &in->size, &in->alloc);
        PROGRESS_MSG("Read input");

        rc = pthread_mutex_unlock(&tq->read_mtx);
        CHECK_THREAD(rc);

        /* no need to free any resources here */
        if (! in->size)
            pthread_exit(NULL);

        /* search for EMPTY output buffer, append, and set status
           to LOADING */
        PROGRESS_START();
        while (1)
        {
            rc = pthread_mutex_lock(&tq->pool_mtx); /* ########## POOL LOCK ########## */
            CHECK_THREAD(rc);
            for (p = 0; p != num_pool; ++p)
                if (tq->out_pool[p].status == EMPTY)
                {
                    in->out = &tq->out_pool[p];
                    in->out->status = LOADING;
                    in->out->next = NULL;
                    
                    if (! tq->out_head)
                        tq->out_head = tq->out_tail = in->out;
                    else
                        tq->out_tail->next = in->out, tq->out_tail = in->out;

                    break;
                }
            rc = pthread_mutex_unlock(&tq->pool_mtx); /* ######### POOL UNLOCK ######### */
            CHECK_THREAD(rc);
            if (in->out)
                break;
            sleep(1); 
            /* one second is long enough not to cause too many
               locks/unlocks, but frequent enough compared to the time
               it takes for the main loop */
        }
        PROGRESS_MSG("Buffer search");


        PROGRESS_START();
        /* load the output buffer. */
        tq->worker(in->worker_par, in->buf, in->size, 
                   &in->out->buf, &in->out->size, &in->out->alloc);
        PROGRESS_MSG("Work");

        mutex_change_status(tq, in->out, FULL);
        in->out = NULL;
        
        while (tq->out_head && tq->out_head->status == FULL)
        {
            PROGRESS_START();
            mutex_change_status(tq, tq->out_head, UNLOADING);
            tq->offload(tq->offload_par, tq->out_head->buf, tq->out_head->size);
            PROGRESS_MSG("Unload");

            rc = pthread_mutex_lock(&tq->pool_mtx);
            CHECK_THREAD(rc);

            tq->out_head->status = EMPTY;
            tq->out_head->size = 0;
            tq->out_head = tq->out_head->next;

            rc = pthread_mutex_unlock(&tq->pool_mtx);
            CHECK_THREAD(rc);
        }
    }
}
