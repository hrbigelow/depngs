#include <gsl/gsl_rng.h>
#include <gsl/gsl_math.h>
#include <sys/timeb.h>
#include <string.h>
#include <ctype.h>
#include <assert.h>

#include "nucleotide_stats.h"
#include "gen_pair_comp.h"
#include "thread_queue.h"

int simp_usage()
{
    fprintf(stderr, 
            "Usage: dep simp n_threads FTYPE input.jpd dist_string max_depth n_rep alpha sim1.pileup sim2.pileup\n"
            "FTYPE is one of 'Sanger (+33)' or 'Solexa (+64)'\n\n"
            "Simulate a pair of pileup files, for each combination of (depth1, depth2, distance, rep)\n"
            "D1 and D2 vary in [0, max_depth], distance is chosen from dist_string, and\n"
            "rep varies in [0, n_rep)\n"
            "alpha gives an alpha value to generate the first in a pair of\n"
            "base composition points.\n"
            "\n"
            "The dist.rdb file has fields of '<position><tab><distance>'\n"
            );
    return 1;
}

struct base_qual {
    char call, qual;
    unsigned strand;
};

/* Generate a base from a normalized composition */
unsigned gen_base(double *comp, gsl_rng *rg)
{
    unsigned i = 0;
    double q = gsl_rng_uniform(rg), cc = comp[0];
    while (q > cc) cc += comp[++i];
    return i;
}

#define MIN(a,b) ((a) < (b) ? (a) : (b))

#define MAX_DIST 100
unsigned qoffset;
float dist[MAX_DIST];
char dist_str[MAX_DIST][20];
unsigned n_dist;
static struct nucleotide_stats stats;

void init_dist(const char *csv)
{
    n_dist = 0;
    float *d = dist;
    char *loc;
    while (*csv != '\0' && n_dist != MAX_DIST)
    {
        if (*csv == ',') ++csv;
        *d = strtof(csv, &loc);
        strncpy(dist_str[n_dist], csv, loc - csv);
        csv = loc;
        ++d;
        ++n_dist;
    }

}


struct base_qual gen_call(unsigned bi, gsl_rng *rg)
{
    size_t code = 0, quality, strand;
    float q = gsl_rng_uniform(rg), cumul = stats.cpd[bi][0];
    while (q > cumul) cumul += stats.cpd[bi][++code];

    struct base_qual bq;
    decode_nucleotide(code, &bq.call, &quality, &strand);
    bq.qual = (char)(quality + qoffset);
    bq.strand = (unsigned)strand;
    return bq;
}

/* position updated by scan */
unsigned long current_end, max_current_end;
unsigned n_loci_per_batch;



struct work_unit {
    unsigned long start, end;
};

struct work_par {
    gsl_rng *rg;
    unsigned long max_depth, n_rep, n_dist;
    double alpha[4];
};

void scan_unit(void *par, unsigned max_bytes)
{
    struct work_unit *w = par;
    w->start = current_end;
    w->end = current_end = MIN(w->start + n_loci_per_batch, max_current_end);
}

/* copy the contents of the par to the buf */
void read_unit(void *par, struct managed_buf *bufs)
{
    struct work_unit *w = par;
    memcpy(bufs[0].buf, par, sizeof(struct work_unit));
    bufs[0].size = w->end != w->start;
}


struct file_pair {
    FILE *fh1, *fh2;
};

void offload(void *par, const struct managed_buf *bufs)
{
    struct file_pair *fp = par;
    fwrite(bufs[0].buf, 1, bufs[0].size, fp->fh1);
    fwrite(bufs[1].buf, 1, bufs[1].size, fp->fh2);
    fflush(fp->fh1);
    fflush(fp->fh2);
}

/* Calculates the inverse of:
   i = 
   (depth1 * wp.max_depth * wp.n_rep * wp.n_dist)
   + (depth2 * wp.n_rep * wp.n_dist)
   + (rep * wp.n_dist)
   + dist_i 
*/
void unpack_unit(struct work_par wp,
                 unsigned long i, 
                 unsigned *dist_i,
                 unsigned *depth1, 
                 unsigned *depth2, 
                 unsigned *rep)
{
    *rep = i % wp.n_rep;
    i /= wp.n_rep;
    *depth2 = i % (wp.max_depth - 1) + 1;
    i /= (wp.max_depth - 1);
    *depth1 = i % (wp.max_depth - 1) + 1;
    i /= (wp.max_depth - 1);
    *dist_i = i % wp.n_dist;
}


void onexit(void *par) { }


/* Generate a set of loci pairs. The general notion is to
   pre-calculate the total number of units of work based on  */
void work(void *par,
          const struct managed_buf *in,
          struct managed_buf *out)
{
    struct work_unit u = *(struct work_unit *)in->buf;
    struct work_par wpar = *(struct work_par *)par;
    unsigned long i;
    unsigned dmax[2], r, refbase_i, call_i, s, b, di;
    struct pair_comp pc;
    struct base_qual bq;
    static char nucs[] = "ACGT";

    for (i = u.start; i != u.end; ++i)
    {
        unpack_unit(wpar, i, &di, &dmax[0], &dmax[1], &r);
        assert(dmax[0] != 0);
        assert(dmax[1] != 0);

        /* generate two random points from the simplex at
           the specified distance. */
        pc = gen_pair_comp(wpar.alpha, dist[di], wpar.rg);
        char calls[] = "ACGTacgt";
        refbase_i = gsl_rng_uniform_int(wpar.rg, 4);
        calls[refbase_i] = '.';
        calls[refbase_i + 4] = ',';

        /* simulate from both, to the specified depth */
        for (s = 0; s != 2; ++s)
        {
            double *c = s == 0 ? pc.c1 : pc.c2;

            ALLOC_GROW(out[s].buf, out[s].size + 50 + dmax[s] + dmax[s], out[s].alloc);
            char *wb = out[s].buf + out[s].size;
            
            strcpy(wb, dist_str[di]);
            wb += strlen(dist_str[di]);
            wb += sprintf(wb, "\t%lu\t", i);
            *wb++ = nucs[refbase_i];
            wb += sprintf(wb, "\t%u\t", dmax[s]);
            wb[dmax[s]] = '\t';

            char *wq = wb + dmax[s] + 1, *we = wq + dmax[s];
            while (wq != we)
            {
                b = gen_base(c, wpar.rg);
                bq = gen_call(b, wpar.rg);
                call_i = b + (bq.strand == NUC_MINUS_STRAND ? 4 : 0);
                *wb++ = calls[call_i];
                *wq++ = bq.qual;
            }
            *wq++ = '\n';
            out[s].size = wq - out[s].buf;
        }
    }
}

int main_simp(int argc, char **argv)
{
    if (argc != 10) return simp_usage();

    unsigned n_threads = (unsigned)atoi(argv[1]);
    char *ftype_string = argv[2];
    char *base_qual_params_file = argv[3];
    char *dist_string = argv[4];
    unsigned max_depth = (unsigned)atoi(argv[5]);
    unsigned n_rep = (unsigned)atoi(argv[6]);
    double av = atof(argv[7]);

    char *sim1_file = argv[8];
    char *sim2_file = argv[9];



    if (strcmp(ftype_string, "Sanger") == 0) qoffset = 33;
    else if (strcmp(ftype_string, "Solexa") == 0) qoffset = 64;
    else fprintf(stderr, "Error: FTYPE must be 'Sanger' or 'Solexa'\n");

    FILE *sim1_fh = fopen(sim1_file, "w");
    if (! sim1_fh)
    {
        fprintf(stderr, "Couldn't open sim1_file %s for writing\n", sim1_file);
        exit(1);
    }

    FILE *sim2_fh = fopen(sim2_file, "w");
    if (! sim2_fh)
    {
        fprintf(stderr, "Couldn't open sim2_file %s for writing\n", sim2_file);
        exit(1);
    }

    nucleotide_stats_initialize(base_qual_params_file, &stats);
    init_dist(dist_string);

    /* timeb millitime; */
    /* ftime(& millitime); */
    /* gsl_rng_set(rg, millitime.millitm); */

    struct work_unit reader_buf[1] = { { 0, 0 } };
    void *reader_par = &reader_buf[0];
    struct work_par *worker_buf = malloc(n_threads * sizeof(struct work_par));
    void **worker_par = malloc(n_threads * sizeof(void *));

    struct file_pair offload_par = { sim1_fh, sim2_fh };
    thread_queue_reader_t reader = { read_unit, scan_unit };

    unsigned t;
    for (t = 0; t != n_threads; ++t)
    {
        worker_buf[t] = (struct work_par){
            gsl_rng_alloc(gsl_rng_taus),
            max_depth,
            n_rep,
            n_dist,
            { av, av, av, av }
        };
        worker_par[t] = &worker_buf[t];
    }

    /* initialize globals */
    current_end = 0;
    max_current_end = (max_depth - 1) * (max_depth - 1) * n_rep * n_dist;
    n_loci_per_batch = 100000;

    struct thread_queue *tqueue =
        thread_queue_init(reader, &reader_par,
                          work, worker_par,
                          offload, &offload_par,
                          onexit,
                          n_threads,
                          n_threads * 10,
                          1,
                          1,
                          2,
                          1e5);

    thread_queue_run(tqueue);
    thread_queue_free(tqueue);

    for (t = 0; t != n_threads; ++t)
        gsl_rng_free(worker_buf[t].rg);

    free(worker_buf);
    free(worker_par);

    fclose(sim1_fh);
    fclose(sim2_fh);

    return 0;
}
