#include "thread_queue.h"

#include <pthread.h>
#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <string.h>

enum buf_status { EMPTY, LOADING, FULL, UNLOADING };

/* managed output buffers.  There will be T + E of these and they will
   be owned by 0 or 1 thread at any given time (but not necessarily
   the same thread each time) */
struct output_node {
    struct output_node *next;
    struct managed_buf *buf;
    enum buf_status status;
};


/* There will be one of these for each thread, and each thread will
   use the same one for its lifespan. */
struct thread_comp_input {
    struct managed_buf *buf;
    struct output_node *out;
    void *scan_info; /* in/out parameter for 'scan'.  also passed to
                        'read' and 'worker' as a means of
                        communicating with them */
    struct thread_queue *tq;
};


/* the main configuration object. there will be just one of these for
   a given program (typically) or for a given task. */
struct thread_queue {
    thread_queue_reader_t reader;
    thread_queue_worker_t *worker;
    thread_queue_offload_t *offload;
    thread_queue_create_t *on_create;
    thread_queue_exit_t *on_exit;

    /* ensures one-at-a-time invocation of reader.scan */
    pthread_mutex_t scan_mtx;

    unsigned n_max_reading;

    /* protects n_reading. */
    pthread_mutex_t read_mtx;

    /* # of threads currently invoking reader.read. */
    unsigned n_reading; 

    /* condition n_reading != n_max_reading */
    pthread_cond_t read_slot_avail;

    void *offload_par;

    /* protects out_head, pool_status,  */
    pthread_mutex_t out_mtx;
    pthread_cond_t out_buf_avail;
    struct output_node *out_pool, *out_head, *out_tail;

    unsigned pool_status[3]; /* number of output nodes in each of the states */


    struct thread_comp_input *input;
    pthread_t *threads;
    unsigned n_threads, n_extra, n_inputs, n_outputs;
};


struct timespec program_start_time;


/* initialize resources */
struct thread_queue *
thread_queue_init(thread_queue_reader_t reader,
                  void **reader_pars,
                  thread_queue_worker_t worker,
                  thread_queue_offload_t offload,
                  void *offload_par,
                  thread_queue_create_t on_create,
                  thread_queue_exit_t on_exit,
                  unsigned n_threads,
                  unsigned n_extra,
                  unsigned n_max_reading,
                  unsigned n_inputs,
                  unsigned n_outputs,
                  unsigned long max_input_mem)
{
    unsigned n_pool = n_threads + n_extra;
    struct thread_queue *tq = malloc(sizeof(struct thread_queue));
    tq->reader = reader;
    tq->n_max_reading = n_max_reading;
    pthread_mutex_init(&tq->read_mtx, NULL);
    tq->n_reading = 0;
    pthread_mutex_init(&tq->scan_mtx, NULL);
    pthread_cond_init(&tq->read_slot_avail, NULL);
    tq->worker = worker;
    tq->offload = offload;
    tq->offload_par = offload_par;
    tq->on_create = on_create;
    tq->on_exit = on_exit;
    tq->out_pool = calloc(n_pool, sizeof(struct output_node));
    tq->out_head = NULL;
    tq->out_tail = NULL;
    memset(tq->pool_status, 0, sizeof(tq->pool_status));
    tq->pool_status[EMPTY] = n_pool;
    tq->input = calloc(n_threads, sizeof(struct thread_comp_input));

    pthread_mutex_init(&tq->out_mtx, NULL);
    pthread_cond_init(&tq->out_buf_avail, NULL);
    tq->threads = malloc(n_threads * sizeof(pthread_t));
    tq->n_threads = n_threads;
    tq->n_extra = n_extra;
    tq->n_inputs = n_inputs;
    tq->n_outputs = n_outputs;

    /* only need initialize the pointer to NULL */
    unsigned p, t, b;
    for (t = 0; t != n_threads; ++t) {
        tq->input[t].buf = malloc(n_inputs * sizeof(struct managed_buf));
        for (b = 0; b != n_inputs; ++b) {
            struct managed_buf *buf = &tq->input[t].buf[b];
            buf->alloc = max_input_mem / n_threads / n_inputs;
            buf->buf = malloc(buf->alloc);
            buf->size = 0;
        }
        tq->input[t].scan_info = reader_pars[t];
        tq->input[t].tq = tq;
    }

    /* by default, use 2% of the memory of the input for the whole
       initial output */
    unsigned out_chunk_size = max_input_mem * 2 / 100 / n_pool + 10;
    for (p = 0; p != n_pool; ++p) {
        tq->out_pool[p].buf = malloc(n_outputs * sizeof(struct managed_buf));
        for (b = 0; b != n_outputs; ++b) {
            tq->out_pool[p].buf[b].alloc = out_chunk_size;
            tq->out_pool[p].buf[b].size = 0;
            tq->out_pool[p].buf[b].buf = malloc(tq->out_pool[p].buf[b].alloc);
        }
    }
    clock_gettime(CLOCK_REALTIME, &program_start_time);

    return tq;
}


/* release resources */
void
thread_queue_free(struct thread_queue *tq)
{
    pthread_mutex_destroy(&tq->out_mtx);
    pthread_mutex_destroy(&tq->scan_mtx);
    pthread_mutex_destroy(&tq->read_mtx);
    pthread_cond_destroy(&tq->read_slot_avail);
    pthread_cond_destroy(&tq->out_buf_avail);

    free(tq->threads);

    unsigned t, b;
    for (t = 0; t != tq->n_threads; ++t) {
        for (b = 0; b != tq->n_inputs; ++b)
            free(tq->input[t].buf[b].buf);
        free(tq->input[t].buf);
    }
    free(tq->input);
    unsigned n_pool = tq->n_threads + tq->n_extra;
    for (t = 0; t != n_pool; ++t) {
        for (b = 0; b != tq->n_outputs; ++b)
            free(tq->out_pool[t].buf[b].buf);
        free(tq->out_pool[t].buf);
    }
    free(tq->out_pool);
    free(tq);
}


static void *
worker_func(void *args);


#define CHECK_THREAD(rc)                                                \
    if (rc) {                                                           \
        fprintf(stderr,                                                 \
                "Thread-related error at %s:%u with return code %i\n",  \
                __FILE__, __LINE__, rc);                                \
        exit(1);                                                        \
    }


/* start things running. */
int
thread_queue_run(struct thread_queue *tq)
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


/* safely change the status flag of the 'out' node. */
void
set_outnode_status(struct thread_queue *tq,
                   struct output_node *out,
                   enum buf_status status)
{
    --tq->pool_status[out->status];
    out->status = status;
    ++tq->pool_status[out->status];
    unsigned b;
    if (status == EMPTY)
        for (b = 0; b != tq->n_outputs; ++b)
            out->buf[b].size = 0;
}


#define ELAPSED_MS \
    ((((end_time).tv_sec * 1000000000 + (end_time).tv_nsec) -           \
      ((beg_time).tv_sec * 1000000000 + (beg_time).tv_nsec)) / 1000000)

#ifdef TQ_DBG
#define PROGRESS_DECL() struct timespec beg_time, end_time;

#define TIME_MS(t) \
    (long)((t).tv_sec * 1e3 + (t).tv_nsec / 1e6)                        \
    - (long)(program_start_time.tv_sec * 1e3                            \
             + program_start_time.tv_nsec / 1e6)

#define PROGRESS_START(category)                                        \
    do {                                                                \
    clock_gettime(CLOCK_REALTIME, &beg_time);                           \
    fprintf(stderr,                                                     \
            "START\t%s\t%li\t%li\t%li\t%li\t%u\t%u\t%u\n",              \
            (category),                                                 \
            TIME_MS(beg_time),                                          \
            TIME_MS(beg_time),                                          \
            par - tq->input,                                            \
            0l,                                                         \
            tq->pool_status[EMPTY],                                     \
            tq->pool_status[LOADING],                                   \
            tq->pool_status[FULL]);                                     \
    fflush(stderr);                                                     \
    } while (0)


    
#define PROGRESS_MSG(category)                                          \
    do {                                                                \
        clock_gettime(CLOCK_REALTIME, &end_time);                       \
        fprintf(stderr,                                                 \
                "END\t%s\t%li\t%li\t%li\t%li\t%u\t%u\t%u\n",            \
                (category),                                             \
                TIME_MS(beg_time),                                      \
                TIME_MS(end_time),                                      \
                par - tq->input,                                        \
                ELAPSED_MS,                                             \
                tq->pool_status[EMPTY],                                 \
                tq->pool_status[LOADING],                               \
                tq->pool_status[FULL]);                                 \
        fflush(stderr);                                                 \
    } while (0)
#else
#define PROGRESS_DECL()
#define PROGRESS_START(category)
#define PROGRESS_MSG(category)
#endif


static void
reserve_out_buffer(struct thread_comp_input *par)
{
    struct thread_queue *tq = par->tq;
    unsigned p, n_pool = tq->n_threads + tq->n_extra;

    for (p = 0; p != n_pool; ++p)
        if (tq->out_pool[p].status == EMPTY) {
            par->out = &tq->out_pool[p];
            --tq->pool_status[EMPTY];
            par->out->status = LOADING;
            ++tq->pool_status[LOADING];
            par->out->next = NULL;
                    
            if (! tq->out_head)
                tq->out_head = tq->out_tail = par->out;
            else
                tq->out_tail->next = par->out, tq->out_tail = par->out;

            break;
        }
}


static void
offload_full_buffers(struct thread_comp_input *par)
{
    struct thread_queue *tq = par->tq;
    while (tq->out_head && tq->out_head->status == FULL) {
        tq->pool_status[tq->out_head->status]--;
        tq->out_head->status = UNLOADING;
        tq->pool_status[tq->out_head->status]++;

        /* once tq->out_head->status is set to UNLOADING, program
           logic ensures it is safe to use tq->out_head->buf */
        pthread_mutex_unlock(&tq->out_mtx);
        tq->offload(tq->offload_par, tq->out_head->buf);
        pthread_mutex_lock(&tq->out_mtx);

        tq->pool_status[tq->out_head->status]--;
        tq->out_head->status = EMPTY;
        tq->pool_status[tq->out_head->status]++;
        
        unsigned b;
        for (b = 0; b != tq->n_outputs; ++b)
            tq->out_head->buf[b].size = 0;
        
        tq->out_head = tq->out_head->next;
        pthread_cond_signal(&tq->out_buf_avail);
    }
}


/* this function is run by the thread */
static void *
worker_func(void *args)
{
    struct thread_comp_input *par = args;
    struct thread_queue *tq = par->tq;

    PROGRESS_DECL();
    tq->on_create();

    unsigned more_input = 1;

    while (more_input) {
        PROGRESS_START("SCANWAIT");
        pthread_mutex_lock(&tq->scan_mtx);
        PROGRESS_MSG("SCANWAIT");

        /* reserve next input range */
        PROGRESS_START("SCAN");
        tq->reader.scan(par->scan_info, par->buf[0].alloc);
        PROGRESS_MSG("SCAN");

        /* reserve the first empty out buffer.  we must do this within
           the scan_mtx to ensure output order.  */
        if (tq->n_outputs) {
            PROGRESS_START("RESERVEOUT");
            pthread_mutex_lock(&tq->out_mtx);
            while (! tq->pool_status[EMPTY])
                pthread_cond_wait(&tq->out_buf_avail, &tq->out_mtx);
            
            reserve_out_buffer(par);
            pthread_mutex_unlock(&tq->out_mtx);
            PROGRESS_MSG("RESERVEOUT");
        }

        pthread_mutex_unlock(&tq->scan_mtx);

        /* wait and reserve reader slot */
        PROGRESS_START("READWAIT");
        pthread_mutex_lock(&tq->read_mtx);
        while (tq->n_reading == tq->n_max_reading)
            pthread_cond_wait(&tq->read_slot_avail, &tq->read_mtx);
        ++tq->n_reading;
        pthread_mutex_unlock(&tq->read_mtx);
        PROGRESS_MSG("READWAIT");

        PROGRESS_START("READ");
        more_input = tq->reader.read(par->scan_info, par->buf);
        PROGRESS_MSG("READ");

        pthread_mutex_lock(&tq->read_mtx);
        --tq->n_reading;
        pthread_mutex_unlock(&tq->read_mtx);
        pthread_cond_signal(&tq->read_slot_avail);

        /* load the output buffer. */
        PROGRESS_START("WORK");
        struct managed_buf *reserved_out_buf = tq->n_outputs ? par->out->buf : NULL;
        tq->worker(par->buf, more_input, par->scan_info, reserved_out_buf);
        PROGRESS_MSG("WORK");

        if (tq->n_outputs) {
            pthread_mutex_lock(&tq->out_mtx);
            set_outnode_status(tq, par->out, FULL);
            par->out = NULL;
            offload_full_buffers(par);
            pthread_mutex_unlock(&tq->out_mtx);
        }
    }
    tq->on_exit();
    pthread_exit(NULL);
}
