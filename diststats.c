#include "metropolis_sampling.h"
#include "dirichlet_diff_cache.h"
#include "dirichlet_points_gen.h"
#include "cache.h"
#include "thread_queue.h"
#include "yepLibrary.h"

#include <stdio.h>
#include <string.h>
#include <gsl/gsl_errno.h>
#include <assert.h>


int diststats_usage()
{
    fprintf(stderr,
            "\nUsage: dep diststats [options] output_stats.rdb\n"
            "Options:\n\n"
            "-1 INT      MAX1, maximum first dirichlet alpha component [100]\n"
            "-2 INT      MAX2, maximum second dirichlet alpha component [10]\n"
            "-y FLOAT    MIN_DIST, min mutation distance (0-1 scale) to call as changed [0.2]\n"
            "-X FLOAT    POST_CONF, confidence for -y [0.99]\n"
            "-Z FLOAT    BETA_CONF, confidence for binomial estimation [0.9999]\n"
            "-f INT      MAX_POINTS, max # of sample point pairs for binomial test [10000]\n"
            "-t INT      number of threads to use [1]\n"
            "-p FLOAT    PRIOR_ALPHA, alpha value for dirichlet prior [0.1]\n"
            "\n"
            "Computes a statistic over a space of pairs of 4D dirichlet\n"
            "distributions.  The first takes alpha values (p+a1, p+a2, p,\n"
            "p), (p = PRIOR_ALPHA).  The second takes alpha values (p+b1,\n"
            "p+b2, p, p).  All combinations of values of a1, a2, b1, b2 are\n"
            "computed, with a1 and b1 in [0, MAX1) a2 and b2 in [0, MAX2), a1\n"
            ">= a2 and b1 >= b2.  For each of these distribution pairs, up to\n"
            "MAX_POINTS pairs of points are drawn, one point at a time from\n"
            "each distribution the euclidean distance is computed for each\n"
            "pair, and marked as a bernoulli success (S) if the distance is\n"
            "less than MIN_DIST, failure (F) otherwise.  The Beta distribution\n"
            "Beta(S + 1/2, F + 1/2) (Jeffrey's test) is used with BETA_CONF\n"
            "and 1-BETA_CONF quantiles as lower and upper bound estimates on\n"
            "the probability of success (BETA_LB, BETA_UB).\n"
            "\n"
            "The final statistic is called as follows:\n"
            "BETA_UB < (1-POST_CONF): UNCHANGED\n"
            "BETA_LB > POST_CONF    : CHANGED\n"
            "otherwise              : AMBIGUOUS\n"
            "\n"
            "This table is then used in dep dist for rapid classification of\n"
            "pairs of loci.\n"
            );
    return 1;
}


struct alpha_limits {
    unsigned max1, max2;
};

struct triplet {
    unsigned a2, b1, b2;
};

/* reader function.  simply return the next quartet of values, encoded
   in the buf. signal end of input by settings buf[0].size = 0 */
void diststats_reader(void *par, struct managed_buf *buf)
{
    ALLOC_GROW(buf[0].buf, sizeof(struct triplet), buf[0].alloc);
    static struct triplet alpha = { 0, 0, 0 };
    struct alpha_limits *lim = par;
    if (alpha.a2 < lim->max2)
    {
        memcpy(buf[0].buf, &alpha, sizeof(alpha));
        buf[0].size = sizeof(struct triplet);
        ++alpha.b2;
        if (alpha.b2 == lim->max2) alpha.b2 = 0, ++alpha.b1;
        if (alpha.b1 == lim->max1) alpha.b1 = 0, ++alpha.a2;
    }
    else buf[0].size = 0;
}


void diststats_worker(void *par,
                      const struct managed_buf *in_bufs,
                      struct managed_buf *out_bufs)
{
    struct binomial_est_params *bep = par;
    struct binomial_est_bounds beb;
    struct triplet t;

    ALLOC_GROW(out_bufs[0].buf, 
               sizeof(beb) + sizeof(t),
               out_bufs[0].alloc);

    memcpy(&t, in_bufs[0].buf, sizeof(struct triplet));
    initialize_est_bounds(t.a2, t.b1, t.b2, bep, &beb);
    memcpy(out_bufs[0].buf, &t, sizeof(t));
    memcpy(out_bufs[0].buf + sizeof(t), &beb, sizeof(beb));
    out_bufs[0].size = sizeof(beb);
}


void diststats_offload(void *par,
                       const struct managed_buf *bufs)
{
    FILE *out_fh = par;
    struct binomial_est_bounds beb;
    struct triplet t;
    memcpy(&t, bufs[0].buf, sizeof(t));
    memcpy(&beb, bufs[0].buf + sizeof(t), sizeof(beb));
    fprintf(out_fh,
            "%i\t%i\t%i\t%i\t%i\t%i\t%i\n", 
            t.a2, t.b1, t.b2,
            beb.ambiguous[0], beb.unchanged[0], beb.unchanged[1], beb.ambiguous[1]);
}



int main_diststats(int argc, char **argv)
{
    size_t n_threads = 1;
    struct alpha_limits reader_par = { 100, 10 };

    float prior_alpha = 0.1;
    struct posterior_settings pset = {
        { prior_alpha, prior_alpha, prior_alpha, prior_alpha },
        10000,
        0.2,
        0.99,
        0.9999
    };

    char c;
    while ((c = getopt(argc, argv, "1:2:t:f:y:X:Z:p:")) >= 0)
    {
        switch(c)
        {
        case '1': reader_par.max1 = (size_t)atof(optarg); break;
        case '2': reader_par.max2 = (size_t)atof(optarg); break;
        case 't': n_threads = (size_t)atof(optarg); break;
        case 'f': pset.max_sample_points = (size_t)atof(optarg); break;
        case 'y': pset.min_dist = atof(optarg); break;
        case 'X': pset.post_confidence = atof(optarg); break;
        case 'Z': pset.beta_confidence = atof(optarg); break;
        case 'p': pset.prior_alpha[0] = 
                pset.prior_alpha[1] = 
                pset.prior_alpha[2] = 
                pset.prior_alpha[3] = atof(optarg); break;
        default: return diststats_usage(); break;
        }
    }
    if (argc - optind != 1) return diststats_usage();

    char *out_file = argv[optind];
    FILE *out_fh = fopen(out_file, "w");
    if (! out_fh)
    {
        fprintf(stderr, "Error: Couldn't open output file %s for writing\n",
                out_file);
        exit(1);
    }

    gsl_set_error_handler_off();

    init_beta(pset.beta_confidence);
    dirichlet_diff_init();

    struct binomial_est_params worker_par;
    worker_par.pset = &pset;
    worker_par.dist[0] = malloc(sizeof(struct distrib_points));
    worker_par.dist[1] = malloc(sizeof(struct distrib_points));
    alloc_distrib_points(worker_par.dist[0], pset.max_sample_points);
    alloc_distrib_points(worker_par.dist[1], pset.max_sample_points);
    worker_par.batch_size = GEN_POINTS_BATCH;

    void *worker_par_ary[] = { &worker_par };
    
    size_t n_extra = n_threads;
    struct thread_queue *tqueue =
        thread_queue_init(diststats_reader, &reader_par,
                          diststats_worker, worker_par_ary,
                          diststats_offload, out_fh,
                          n_threads,
                          n_extra,
                          1,
                          1,
                          1e6);

    enum YepStatus status = yepLibrary_Init();
    assert(status == YepStatusOk);

    thread_queue_run(tqueue);
    thread_queue_free(tqueue);
    dirichlet_diff_free();
    free_distrib_points(worker_par.dist[0]);
    free_distrib_points(worker_par.dist[1]);
    free(worker_par.dist[0]);
    free(worker_par.dist[1]);
    return 0;
}
