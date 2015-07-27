#include "thread_queue.h"

#include <pthread.h>
#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <string.h>

enum buf_status {
    EMPTY,
    LOADING,
    FULL
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
    thread_queue_create_t *on_create;
    thread_queue_exit_t *on_exit;

    pthread_mutex_t io_mtx;
    pthread_cond_t reader_free_cond;

    /* protected by io_mtx */
    void **reader_par;
    unsigned *reader_in_use, n_readers, n_readers_free;

    void *offload_par;

    pthread_mutex_t out_mtx;
    pthread_cond_t out_free_cond;
    struct output_node *out_pool, *out_head, *out_tail;

    unsigned pool_status[3]; /* number of output nodes in each of the states */
    struct thread_comp_input *input;
    pthread_t *threads;
    unsigned n_threads, n_extra, n_inputs, n_outputs;
};


struct timespec program_start_time;


/* initialize resources */
struct thread_queue *
thread_queue_init(thread_queue_reader_t reader, void **reader_par,
                  thread_queue_worker_t worker,
                  thread_queue_offload_t offload, void *offload_par,
                  thread_queue_create_t on_create,
                  thread_queue_exit_t on_exit,
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
    pthread_mutex_init(&tq->io_mtx, NULL);
    pthread_cond_init(&tq->reader_free_cond, NULL);
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
    pthread_cond_init(&tq->out_free_cond, NULL);
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
            tq->input[t].buf[b].alloc = max_input_mem / n_threads / n_inputs;
            tq->input[t].buf[b].buf = malloc(tq->input[t].buf[b].alloc);
            tq->input[t].buf[b].size = 0;
        }
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


/* release resources */
void
thread_queue_free(struct thread_queue *tq)
{
    pthread_mutex_destroy(&tq->out_mtx);
    pthread_mutex_destroy(&tq->io_mtx);
    pthread_cond_destroy(&tq->reader_free_cond);
    pthread_cond_destroy(&tq->out_free_cond);

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
}


/* safely change the status flag of the 'out' node. */
void
set_outnode_status(struct thread_queue *tq,
                   struct output_node *out,
                   enum buf_status status)
{
    pthread_mutex_lock(&tq->out_mtx);
    --tq->pool_status[out->status];
    out->status = status;
    ++tq->pool_status[out->status];
    unsigned b;
    if (status == EMPTY)
        for (b = 0; b != tq->n_outputs; ++b)
            out->buf[b].size = 0;

    pthread_mutex_unlock(&tq->out_mtx);
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
            in - tq->input,                                             \
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
                in - tq->input,                                         \
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


/* this function is run by the thread */
static void *
worker_func(void *args)
{
    struct thread_comp_input *in = args;
    struct thread_queue *tq = in->tq;
    unsigned p, n_pool = tq->n_threads + tq->n_extra;

    PROGRESS_DECL();
    tq->on_create();

    while (1) {
        PROGRESS_START("WAIT");
        pthread_mutex_lock(&tq->io_mtx);
        while (! tq->n_readers_free)
            pthread_cond_wait(&tq->reader_free_cond, &tq->io_mtx);

        /* reserve first free reader */
        unsigned r;
        for (r = 0; r != tq->n_readers; ++r)
            if (! tq->reader_in_use[r]) break;
        
        tq->reader_in_use[r] = 1;
        --tq->n_readers_free;
        PROGRESS_MSG("WAIT");

        /* reserve next input range */
        PROGRESS_START("SCAN");
        tq->reader.scan(tq->reader_par[r], in->buf[0].alloc);
        PROGRESS_MSG("SCAN");

        /* reserve the first empty out buffer */
        pthread_mutex_lock(&tq->out_mtx);
        if (! tq->pool_status[EMPTY])
            pthread_cond_wait(&tq->out_free_cond, &tq->out_mtx);

        for (p = 0; p != n_pool; ++p)
            if (tq->out_pool[p].status == EMPTY) {
                in->out = &tq->out_pool[p];
                --tq->pool_status[EMPTY];
                in->out->status = LOADING;
                ++tq->pool_status[LOADING];
                in->out->next = NULL;
                    
                if (! tq->out_head)
                    tq->out_head = tq->out_tail = in->out;
                else
                    tq->out_tail->next = in->out, tq->out_tail = in->out;

                break;
            }
        pthread_mutex_unlock(&tq->out_mtx);
        pthread_mutex_unlock(&tq->io_mtx);

        PROGRESS_START("READ");
        tq->reader.read(tq->reader_par[r], in->buf);
        PROGRESS_MSG("READ");

        pthread_mutex_lock(&tq->io_mtx);
        tq->reader_in_use[r] = 0;
        ++tq->n_readers_free;
        pthread_cond_signal(&tq->reader_free_cond);
        pthread_mutex_unlock(&tq->io_mtx);

        /* no need to free any resources here */
        unsigned sz = 0, b;
        for (b = 0; b != tq->n_inputs; ++b)
            if (in->buf[b].size) {
                sz = 1;
                break;
            }
        if (sz == 0) {
            /* need to free up output buffer. setting it to full just
               has the effect of holding onto one extra buffer
               unnecessarily. when it comes time to output it, the
               length will be zero, so it will have no effect. */
            set_outnode_status(tq, in->out, FULL);
            tq->on_exit();
            pthread_exit(NULL);
        }

        /* load the output buffer. */
        PROGRESS_START("WORK");
        tq->worker(in->buf, in->out->buf);
        PROGRESS_MSG("WORK");

        set_outnode_status(tq, in->out, FULL);
        in->out = NULL;
        
        pthread_mutex_lock(&tq->out_mtx);
        while (tq->out_head && tq->out_head->status == FULL) {
            tq->offload(tq->offload_par, tq->out_head->buf);
            --tq->pool_status[tq->out_head->status];
            tq->out_head->status = EMPTY;
            ++tq->pool_status[tq->out_head->status];
            
            for (b = 0; b != tq->n_outputs; ++b)
                tq->out_head->buf[b].size = 0;

            tq->out_head = tq->out_head->next;
        }
        pthread_mutex_unlock(&tq->out_mtx);
    }
}
