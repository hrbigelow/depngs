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
#ifndef _THREAD_QUEUE_H
#define _THREAD_QUEUE_H

#include <unistd.h>
#include "cache.h"

/* the main configuration object. opaque type using pimpl pattern */
struct thread_queue;

/* this client-provided reader function is expected to manage its
   buf/size/alloc space as needed.  The thread_queue mechanism will
   pass NULL, 0, 0 as arguments for the first call to reader. The
   thread_queue library will de-allocate once work is done
   however.  */
typedef void (thread_queue_reader_t)(void *par, struct managed_buf *bufs);

/* the client-provided worker consumes the input produced from the
   reader and outputs it into out_buf, managing out_size and
   out_alloc.*/
typedef void (thread_queue_worker_t)(void *par, 
                                     const struct managed_buf *in_bufs,
                                     struct managed_buf *out_bufs);

/* this function will be called on each output chunk in input order,
   and as soon as the output chunk is finished writing. buf and size
   are the buffer and size of output produced by the writer function.
   par controls the behavior of this offloading function.  once the
   offload is called, the buffer can be re-used by a new thread.
   */
typedef void (thread_queue_offload_t)(void *par, 
                                      const struct managed_buf *bufs);

/* initialize resources */
struct thread_queue *
thread_queue_init(thread_queue_reader_t reader,
                  void *reader_par,
                  thread_queue_worker_t worker,
                  void *worker_par,
                  thread_queue_offload_t offload,
                  void *offload_par,
                  size_t num_threads,
                  size_t num_extra_in_pool,
                  size_t num_inputs,
                  size_t num_outputs,
                  size_t max_input_mem);


/* start things running. */
int thread_queue_run(struct thread_queue *tq);


/* release resources */
void thread_queue_free(struct thread_queue *tq);

#endif /* THREAD_QUEUE_H */
