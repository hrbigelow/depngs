#include "metropolis_sampling.h"
#include "dir_diff_cache.h"
#include "dir_points_gen.h"
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
            "-s INT      RANGE_START, b1 start value [0]\n"
            "-e INT      RANGE_END, b1 end value [MAX1]\n"
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
            "\n"
            "To help parallelize, the full range of b1 values [0, MAX1] may be\n"
            "broken up into [RANGE_START, RANGE_END) half-open sub-intervals\n"
            "and then merged.\n"
            );
    return 1;
}


struct alpha_limits {
    unsigned max1, max2;
    unsigned b1_beg, b1_end;
    struct alpha_pair alpha;
};

/* reader function.  simply return the next triplet of values, encoded
   in the buf. signal end of input by settings buf[0].size = 0
   space over a2:[0, max2), b1:[0, max1), b2:[0, max2)
 */
void diststats_reader(void *par, struct managed_buf *buf)
{
    ALLOC_GROW(buf[0].buf, sizeof(struct alpha_pair), buf[0].alloc);
    struct alpha_limits *lim = par;
    if (lim->alpha.b1 != lim->b1_end)
    {
        memcpy(buf[0].buf, &lim->alpha, sizeof(lim->alpha));
        buf[0].size = sizeof(struct alpha_pair);
        ++lim->alpha.b2;
        if (lim->alpha.b2 == lim->max2) lim->alpha.b2 = 0, ++lim->alpha.b1;
    }
    else buf[0].size = 0;
}


struct work_unit {
    unsigned a2;
    struct alpha_pair p;
    struct binomial_est_bounds beb;
};

void diststats_worker(void *par,
                      const struct managed_buf *in_bufs,
                      struct managed_buf *out_bufs)
{
    struct binomial_est_params *bep = par;
    struct binomial_est_bounds beb;
    struct alpha_pair bpair;

    size_t space = bep->max2 * sizeof(struct work_unit);
    ALLOC_GROW(out_bufs[0].buf, space, out_bufs[0].alloc);

    memcpy(&bpair, in_bufs[0].buf, sizeof(struct alpha_pair));
    char *out = out_bufs[0].buf;
    struct work_unit unit;

    unsigned a2;
    for (a2 = 0; a2 != bep->max2; ++a2)
    {
        initialize_est_bounds(a2, bpair.b1, bpair.b2, bep, &beb);
        unit = (struct work_unit){ a2, bpair, beb };
        memcpy(out, &unit, sizeof(unit));
        out += sizeof(unit);
    }
    out_bufs[0].size = space;
    printf("Finished bounds for (b1, b2) =  %u, %u\n", bpair.b1, bpair.b2);
}


struct offload_par {
    FILE *out_fh;
    unsigned max2;
};

void diststats_offload(void *par,
                       const struct managed_buf *bufs)
{
    
    struct offload_par *offpar = par;
    struct work_unit unit;
    unsigned a2;
    char *in = bufs[0].buf;
    for (a2 = 0; a2 != offpar->max2; ++a2)
    {
        memcpy(&unit, in, sizeof(unit));
        in += sizeof(unit);
        write_diststats_line(offpar->out_fh, a2, 
                             unit.p.b1, unit.p.b2, &unit.beb);
    }
}



int main_diststats(int argc, char **argv)
{
    size_t n_threads = 1;
    struct alpha_limits reader_par = { 
        .max1 = 100,
        .max2 = 10, 
        .b1_beg = 0, 
        .b1_end = 100
    };

    float prior_alpha = 0.1;
    struct posterior_settings pset = {
        { prior_alpha, prior_alpha, prior_alpha, prior_alpha },
        10000,
        0.2,
        0.99,
        0.9999
    };

    char c;
    while ((c = getopt(argc, argv, "1:2:t:f:y:X:Z:p:s:e:")) >= 0)
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
        case 's': reader_par.b1_beg = (size_t)atof(optarg); break;
        case 'e': reader_par.b1_end = (size_t)atof(optarg); break;
        default: return diststats_usage(); break;
        }
    }
    if (argc - optind != 1) return diststats_usage();

    /* initialize the counter with the proper starting value */
    reader_par.alpha = (struct alpha_pair){ .b1 = reader_par.b1_beg, .b2 = 0 };

    char *out_file = argv[optind];
    FILE *out_fh = fopen(out_file, "w");
    if (! out_fh)
    {
        fprintf(stderr, "Error: Couldn't open output file %s for writing\n",
                out_file);
        exit(1);
    }

    gsl_set_error_handler_off();

    setvbuf(stdout, NULL, _IONBF, 0);
    printf("\n"); /* So progress messages don't interfere with shell prompt. */
    
    printf("Precomputing confidence interval statistics...");
    init_beta(pset.beta_confidence, n_threads);
    printf("done.\n");

    dirichlet_diff_init(reader_par.max1, reader_par.max2);

    pset.max_sample_points += GEN_POINTS_BATCH - (pset.max_sample_points % GEN_POINTS_BATCH);

    struct binomial_est_params *worker_par = malloc(n_threads * sizeof(struct binomial_est_params));
    void **worker_par_ptrs = malloc(n_threads * sizeof(void *));

    struct offload_par out_par = { out_fh, reader_par.max2 };

    unsigned t;
    for (t = 0; t != n_threads; ++t)
    {
        worker_par_ptrs[t] = &worker_par[t];
        worker_par[t] = (struct binomial_est_params){ 
            .pset = &pset,
            .dist = { malloc(sizeof(struct distrib_points)),
                      malloc(sizeof(struct distrib_points)) },
            .max1 = reader_par.max1,
            .max2 = reader_par.max2,
            .batch_size = GEN_POINTS_BATCH
        };
        alloc_distrib_points(worker_par[t].dist[0], pset.max_sample_points);
        alloc_distrib_points(worker_par[t].dist[1], pset.max_sample_points);
    }

    write_diststats_header(out_fh, pset, reader_par.max1, reader_par.max2);

    size_t n_extra = n_threads;
    struct thread_queue *tqueue =
        thread_queue_init(diststats_reader, &reader_par,
                          diststats_worker, worker_par_ptrs,
                          diststats_offload, &out_par,
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
    for (t = 0; t != n_threads; ++t)
    {
        free_distrib_points(worker_par[t].dist[0]);
        free_distrib_points(worker_par[t].dist[1]);
    }
    free(worker_par);
    free(worker_par_ptrs);
    return 0;
}
