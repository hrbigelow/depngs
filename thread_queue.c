#include "thread_queue.h"

#include <pthread.h>
#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <string.h>

enum buf_status {
    EMPTY,
    LOADING,
    FULL,
    UNLOADING
};

/* managed output buffers.  There will be T + E of these and they will
   be owned by 0 or 1 thread at any given time (but not necessarily
   the same thread each time) */
struct output_node {
    struct output_node *next;
    struct managed_buf *buf;
    size_t n_buf; /* unused and uninitialized. */
    enum buf_status status;
};


/* There will be one of these for each thread, and each thread will
   use the same one for its lifespan. */
struct thread_comp_input {
    void *worker_par; /* generic parameters to be passed to the worker */
    struct managed_buf *buf;
    size_t n_buf;
    struct output_node *out;
    struct thread_queue *tq;
};


/* the main configuration object. there will be just one of these for
   a given program (typically) or for a given task. */
struct thread_queue {
    thread_queue_reader_t reader;
    thread_queue_worker_t *worker;
    thread_queue_offload_t *offload;
    thread_queue_exit_t *onexit;

    /* parameters to be used by the reader within the read_mtx lock */
    void **reader_par;
    void *global_reader_state; /* allows threads to communicate reader
                                  progress. */
    unsigned *reader_in_use, n_readers, n_readers_free;
    pthread_mutex_t read_mtx;
    pthread_cond_t read_cond;
    void *offload_par;
    struct output_node *out_pool, *out_head, *out_tail;
    unsigned pool_status[4]; /* number of output nodes in each of the states */
    struct thread_comp_input *input;
    pthread_mutex_t pool_mtx;
    pthread_t *threads;
    unsigned n_threads, n_extra, n_inputs, n_outputs;
};


struct timespec program_start_time;


/* initialize resources */
struct thread_queue *
thread_queue_init(thread_queue_reader_t reader, void **reader_par,
                  thread_queue_worker_t worker, void **worker_par,
                  thread_queue_offload_t offload, void *offload_par,
                  thread_queue_exit_t onexit,
                  void *global_reader_state,
                  unsigned n_threads,
                  unsigned n_extra,
                  unsigned n_readers,
                  unsigned n_inputs,
                  unsigned n_outputs,
                  unsigned long max_input_mem)
{
    unsigned n_pool = n_threads + n_extra;
    struct thread_queue *tq = malloc(sizeof(struct thread_queue));
    tq->reader = reader;
    tq->reader_par = reader_par;
    tq->n_readers = n_readers;
    tq->n_readers_free = n_readers;
    tq->reader_in_use = calloc(n_readers, sizeof(tq->reader_in_use[0]));
    pthread_mutex_init(&tq->read_mtx, NULL);
    pthread_cond_init(&tq->read_cond, NULL);
    tq->worker = worker;
    tq->offload = offload;
    tq->offload_par = offload_par;
    tq->onexit = onexit;
    tq->out_pool = calloc(n_pool, sizeof(struct output_node));
    tq->out_head = NULL;
    tq->out_tail = NULL;
    memset(tq->pool_status, 0, sizeof(tq->pool_status));
    tq->pool_status[EMPTY] = n_pool;
    tq->global_reader_state = global_reader_state;
    tq->input = calloc(n_threads, sizeof(struct thread_comp_input));

    pthread_mutex_init(&tq->pool_mtx, NULL);
    pthread_mutex_init(&tq->read_mtx, NULL);
    tq->threads = malloc(n_threads * sizeof(pthread_t));
    tq->n_threads = n_threads;
    tq->n_extra = n_extra;
    tq->n_inputs = n_inputs;
    tq->n_outputs = n_outputs;

    /* only need initialize the pointer to NULL */
    unsigned p, t, b;
    for (t = 0; t != n_threads; ++t)
    {
        tq->input[t].worker_par = worker_par[t];
        tq->input[t].buf = malloc(n_inputs * sizeof(struct managed_buf));
        for (b = 0; b != n_inputs; ++b)
        {
            tq->input[t].buf[b].alloc = max_input_mem / n_threads / n_inputs;
            tq->input[t].buf[b].buf = malloc(tq->input[t].buf[b].alloc);
            tq->input[t].buf[b].size = 0;
        }
        tq->input[t].tq = tq;
    }

    /* by default, use 2% of the memory of the input for the whole
       initial output */
    unsigned out_chunk_size = max_input_mem * 2 / 100 / n_pool + 10;
    for (p = 0; p != n_pool; ++p)
    {
        tq->out_pool[p].buf = malloc(n_outputs * sizeof(struct managed_buf));
        for (b = 0; b != n_outputs; ++b)
        {
            tq->out_pool[p].buf[b].alloc = out_chunk_size;
            tq->out_pool[p].buf[b].size = 0;
            tq->out_pool[p].buf[b].buf = malloc(tq->out_pool[p].buf[b].alloc);
        }
    }
    clock_gettime(CLOCK_REALTIME, &program_start_time);

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

    unsigned t;
    int rc;
    for (t = 0; t != tq->n_threads; ++t) {
        rc = pthread_create(&tq->threads[t], NULL, worker_func, &tq->input[t]);
        CHECK_THREAD(rc);
    }

    for (t = 0; t != tq->n_threads; ++t) {
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

    unsigned t, b;
    for (t = 0; t != tq->n_threads; ++t)
    {
        for (b = 0; b != tq->n_inputs; ++b)
            free(tq->input[t].buf[b].buf);
        free(tq->input[t].buf);
    }
    free(tq->input);
    unsigned n_pool = tq->n_threads + tq->n_extra;
    for (t = 0; t != n_pool; ++t)
    {
        for (b = 0; b != tq->n_outputs; ++b)
            free(tq->out_pool[t].buf[b].buf);
        free(tq->out_pool[t].buf);
    }
    free(tq->out_pool);
}


/* safely change the status flag of the 'out' node. */
void mutex_change_status(struct thread_queue *tq,
                         struct output_node *out,
                         enum buf_status status)
{
    int rc = pthread_mutex_lock(&tq->pool_mtx);
    CHECK_THREAD(rc);
    --tq->pool_status[out->status];
    out->status = status;
    ++tq->pool_status[out->status];
    unsigned b;
    if (status == EMPTY)
        for (b = 0; b != tq->n_outputs; ++b)
            out->buf[b].size = 0;

    rc = pthread_mutex_unlock(&tq->pool_mtx);
    CHECK_THREAD(rc);
}


#define ELAPSED_MS \
    ((((end_time).tv_sec * 1000000000 + (end_time).tv_nsec) -           \
      ((beg_time).tv_sec * 1000000000 + (beg_time).tv_nsec)) / 1000000)

#ifdef _THREAD_QUEUE_DEBUG
#define PROGRESS_DECL() struct timespec beg_time, end_time;

#define TIME_MS(t) \
    (long)((t).tv_sec * 1e3 + (t).tv_nsec / 1e6)                        \
    - (long)(program_start_time.tv_sec * 1e3                            \
             + program_start_time.tv_nsec / 1e6)

#define PROGRESS_START(category)                                        \
    do {                                                                \
    clock_gettime(CLOCK_REALTIME, &beg_time);                           \
    fprintf(stderr,                                                     \
            "START\t%s\t%li\t%li\t%li\t%li\t%u\t%u\t%u\t%u\n",          \
            (category),                                                 \
            TIME_MS(beg_time),                                          \
            TIME_MS(beg_time),                                          \
            in - tq->input,                                             \
            0l,                                                         \
            tq->pool_status[EMPTY],                                     \
            tq->pool_status[LOADING],                                   \
            tq->pool_status[FULL],                                      \
            tq->pool_status[UNLOADING]);                                \
    fflush(stderr);                                                     \
    } while (0)


    
#define PROGRESS_MSG(category)                                          \
    do {                                                                \
        clock_gettime(CLOCK_REALTIME, &end_time);                       \
        fprintf(stderr,                                                 \
                "END\t%s\t%li\t%li\t%li\t%li\t%u\t%u\t%u\t%u\n",        \
                (category),                                             \
                TIME_MS(beg_time),                                      \
                TIME_MS(end_time),                                      \
                in - tq->input,                                         \
                ELAPSED_MS,                                             \
                tq->pool_status[EMPTY],                                 \
                tq->pool_status[LOADING],                               \
                tq->pool_status[FULL],                                  \
                tq->pool_status[UNLOADING]);                            \
        fflush(stderr);                                                 \
    } while (0)
#else
#define PROGRESS_DECL()
#define PROGRESS_START(category)
#define PROGRESS_MSG(category)
#endif


/* this function is run by the thread */
static void *worker_func(void *args)
{
    struct thread_comp_input *in = args;
    struct thread_queue *tq = in->tq;
    unsigned p, n_pool = tq->n_threads + tq->n_extra;

    PROGRESS_DECL();

    while (1)
    {
        int rc;
        /* read chunk into input buffer.  this step can be done
           independent of any locking of the out-buffer pool. */
        PROGRESS_START("WAIT");
        (void)pthread_mutex_lock(&tq->read_mtx);

        /* wait for a reader to free up. */
        while (! tq->n_readers_free)
            pthread_cond_wait(&tq->read_cond, &tq->read_mtx);

        /* find the first free reader */
        unsigned r;
        for (r = 0; r != tq->n_readers; ++r)
            if (! tq->reader_in_use[r]) break;
        
        tq->reader_in_use[r] = 1;
        --tq->n_readers_free;

        tq->reader.get_global_state(tq->reader_par[r], tq->global_reader_state);
        PROGRESS_MSG("WAIT");

        PROGRESS_START("SCAN");
        tq->reader.scan(tq->reader_par[r], in->buf[0].alloc);
        PROGRESS_MSG("SCAN");

        tq->reader.set_global_state(tq->reader_par[r], tq->global_reader_state);
        (void)pthread_mutex_unlock(&tq->read_mtx);

        PROGRESS_START("READ");
        tq->reader.read(tq->reader_par[r], in->buf);
        PROGRESS_MSG("READ");

        (void)pthread_mutex_lock(&tq->read_mtx);
        tq->reader_in_use[r] = 0;
        ++tq->n_readers_free;
        (void)pthread_mutex_unlock(&tq->read_mtx);

        pthread_cond_signal(&tq->read_cond);

        /* no need to free any resources here */
        unsigned sz = 0, b;
        for (b = 0; b != tq->n_inputs; ++b)
            if (in->buf[b].size) 
            {
                sz = 1;
                break;
            }
        if (sz == 0)
        {
            tq->onexit(in->worker_par);
            pthread_exit(NULL);
        }
        else 

        /* search for EMPTY output buf, append, set status to
           LOADING */
        /* PROGRESS_START("BS"); */
        while (1)
        {
            rc = pthread_mutex_lock(&tq->pool_mtx); /* ########## POOL LOCK ########## */
            CHECK_THREAD(rc);
            for (p = 0; p != n_pool; ++p)
                if (tq->out_pool[p].status == EMPTY)
                {
                    in->out = &tq->out_pool[p];
                    --tq->pool_status[in->out->status];
                    in->out->status = LOADING;
                    ++tq->pool_status[in->out->status];
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
        /* PROGRESS_MSG("BS"); */

        /* load the output buffer. */
        PROGRESS_START("WORK");
        tq->worker(in->worker_par, in->buf, in->out->buf);
        PROGRESS_MSG("WORK");

        mutex_change_status(tq, in->out, FULL);
        in->out = NULL;
        
        /* PROGRESS_START("U"); */
        while (tq->out_head && tq->out_head->status == FULL)
        {
            mutex_change_status(tq, tq->out_head, UNLOADING);
            tq->offload(tq->offload_par, tq->out_head->buf);

            rc = pthread_mutex_lock(&tq->pool_mtx);
            CHECK_THREAD(rc);

            --tq->pool_status[tq->out_head->status];
            tq->out_head->status = EMPTY;
            ++tq->pool_status[tq->out_head->status];
            
            for (b = 0; b != tq->n_outputs; ++b)
                tq->out_head->buf[b].size = 0;

            tq->out_head = tq->out_head->next;

            rc = pthread_mutex_unlock(&tq->pool_mtx);
            CHECK_THREAD(rc);
        }
        /* PROGRESS_MSG("U"); */
    }
}
