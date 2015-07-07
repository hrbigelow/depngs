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

/* 'read' uses par to instruct it which ranges of which files to read
   into bufs. thread_queue_free deallocates these buffers.  'scan'
   scans the same files described by par. */
typedef struct {
    void (*read)(void *par, struct managed_buf *bufs);
    void (*scan)(void *par, unsigned max_bytes);
} thread_queue_reader_t;

/* the client-provided worker consumes the input (one or more in_bufs)
   produced from the reader and outputs it into one or more out_bufs,
   managing out_size and out_alloc.  the number of in_bufs and
   out_bufs need not match, and the relation between them is
   determined by the worker logic. */
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


typedef void (thread_queue_create_t)();

typedef void (thread_queue_exit_t)();

/* initialize resources.

   -- void *reader_par must hold the address of an array of n_readers
      addresses, each an address to a struct chosen by the
      user. thread_queue maintains a separate 'in_use' flag for each
      reader.

   -- void *offload_par must hold the address of a single struct,
      which will be cast to that struct by the offload function.

   -- worker_par must hold the address of an array of n_threads
      addresses, each an address to a struct to be used for each
      worker function.  The struct is chosen by the user.
*/
struct thread_queue *
thread_queue_init(thread_queue_reader_t reader, void **reader_par,
                  thread_queue_worker_t worker,
                  thread_queue_offload_t offload, void *offload_par,
                  thread_queue_create_t on_create,
                  thread_queue_exit_t on_exit,
                  unsigned n_threads,
                  unsigned n_extra_in_pool,
                  unsigned n_readers,
                  unsigned n_inputs,
                  unsigned n_outputs,
                  unsigned long max_input_mem);


/* start things running. */
int thread_queue_run(struct thread_queue *tq);


/* release resources */
void thread_queue_free(struct thread_queue *tq);

#endif /* THREAD_QUEUE_H */
