/* 
   Provides a mechanism for multi-threaded processing of a large input
   file.  Each thread reads the next available chunk, processes it and
   produces output in its own buffer.  The mechanism keeps track of
   the order of each chunk and when each is ready for output, as well
   as manages the threads and their resources. 

   Importantly, the N threads have access to an N + E (extra) pool of
   output chunk buffers so that threads that finish their chunk out of
   order may start working on the next available chunk, even before
   the oldest input chunk is ready for output.

   This allows all threads to be working at all times, while still
   producing the output chunks in the same order as the input chunks
   are encountered.

*/

/* the main configuration object. this is intentionally */
struct thread_queue;

/* initialize resources */
struct thread_queue *
thread_queue_init(void *(*worker)(void *),
                  void (*reader)(FILE *fh),
                  size_t num_threads,
                  size_t num_extra_in_pool);


/* start things running. */
int thread_queue_run(struct thread_queue *tq);


/* release resources */
int thread_queue_free(struct thread_queue *tq);
