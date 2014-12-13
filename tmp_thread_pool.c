
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



static int new_query = 1;

/* cast to ptrdiff_t so we can do signed-comparison, even though these
   are always positive. */
#define CHUNK_SIZE() (ptrdiff_t)(base_chunk_size + max_pileup_line_size)
#define FILE_SPAN() (ptrdiff_t)(end_off - start_off)

/* already a ptrdiff_t type */
#define BASE_LEFT() (base_end - write_ptr)


/* */
struct pileup_range_reader {
    FILE *fh;
    struct file_bsearch_index *ix; /* the current index into the file */
    struct locus_range *q, *qend; /* collection of regions desired */
    int new_query;
    char *line_buf;
    size_t max_line_size;
    size_t base_chunk_size;
};


/* should we assume that we are starting out with an allocated buffer?
   secondly, though this may seem strange, this function is passed a
   void *par, which encapsulates all of the state information required
   to read more of the input.  In the simplest case, this will just be
   a FILE *, but in more complex cases, it would also include perhaps
   an index or some other structure needed to specify the input more
   specifically. */

void read_more(void *par, char **in_buf, char **in_end, size_t *in_alloc)
{

    /* allocate a reasonable size for the output buffer */
    size_t t;
    size_t out_buf_size = 1e7;
    off_t start_off, end_off;
    struct pileup_range_reader *rr = par;

    
    /* redo the input strategy assuming off_index input */
    while (rr->q != rr->qend)
    {
        /* kludgy re-use of q != qend test.  is there a better way to
           write this? */
        while (BASE_LEFT() > 0 && rr->q != rr->qend)
        {
            /* find file offsets for current query, creating index
               nodes and updating ix in the process */
            if (rr->new_query)
            {
                rr->ix = find_loose_index(rr->ix, rr->q->beg, rr->fh);
                start_off = off_lower_bound(rr->ix, rr->q->beg);
                rr->ix = find_loose_index(rr->ix, rr->q->end, rr->fh);
                end_off = off_upper_bound(rr->ix, rr->q->end);

                fseeko(rr->fh, start_off, SEEK_SET);
            }

            /* fill the buffer as much as possible with the next query range.
               afterwards, q now points to the next range to retrieve. */
            if (FILE_SPAN() < BASE_LEFT())
            {
                write_ptr += fread(write_ptr, 1, FILE_SPAN(), rr->fh);
                ++rr->q;
                rr->new_query = 1;
            }
            else
            {
                /* partially consume the query range.  read up to the
                   base buffer, then read the next line fragment.
                   realloc both chunk_buf and line_buf as necessary */
                write_ptr += fread(write_ptr, 1, BASE_LEFT(), rr->fh);

                size_t oldmax = max_pileup_line_size;
                ssize_t line_length = getline(&line_buf, &max_pileup_line_size, rr->fh);
                assert(line_length >= 0);

                if (oldmax != max_pileup_line_size)
                {
                    /* should happen very rarely.  if so, realloc and
                       update buffers. */
                    size_t write_pos = write_ptr - chunk_buf;
                    chunk_buf = (char *)realloc(chunk_buf, CHUNK_SIZE());
                    base_end = chunk_buf + base_chunk_size;
                    write_ptr = chunk_buf + write_pos;
                }
                strcpy(write_ptr, line_buf);
                write_ptr += line_length;
                start_off = ftell(pileup_input_fh);

                /* update q->beg to the position just after the last line read */
                char *last_line = (char *)memrchr(chunk_buf, '\n', write_ptr - chunk_buf - 1) + 1;
                rr->q->beg = init_locus(last_line);
                ++rr->q->beg.lo; /* lo is the locus position, hi is the contig */
                rr->new_query = 0;
            }
        }
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



