#ifndef _COMP_WORKER_H
#define _COMP_WORKER_H

#include <vector>
#include <fstream>

#include "pileup_tools.h"
#include "defs.h"

extern "C" {
#include "thread_queue.h"
#include "ordering.h"
#include "metropolis_sampling.h"
#include "cache.h"
#include "binomial_est.h"
}


/* With 100 sample points */
#define MAX_SAMPLE_POINTS 100



/* instantiate one of these for each thread and each sample.  holds
   the information for the locus currently being processed in this
   thread and sample. */
struct locus_sampling
{
    PileupSummary locus;
    bool is_next;
    unsigned char dist_printed;
    /* double proposal_alpha[NUM_NUCS]; */
    struct points_gen pgen;
    struct double_buf points, weights;
    /* size_t autocor_offset; */
    char *current, *end;
    pair_ordering locus_ord;
};

/* attributes intrinsic to one sample */
struct sample_attributes
{
    char label_string[100];
    struct nucleotide_stats nuc_stats;
    FILE *fh;
};


/* instantiate one of these for each thread x sample (but there is
   only one sample being processed) */
struct comp_worker_input
{
    struct sample_attributes sample_atts;
    struct posterior_settings pset;

    double quantiles[MAX_NUM_QUANTILES];
    double quantile_values[MAX_NUM_QUANTILES];
    size_t n_quantiles;

    /* if any non-reference base has its test_quantile greater than
       min_quantile_value it will be reported. */
    float test_quantile, min_test_quantile_value; 
};

/* */
#define GEN_POINTS_BATCH 4

/* the following four functions can be used with binomial_est's struct
   points_gen. */
void gen_dirichlet_points_wrapper(const void *par, double *points);

void calc_post_to_dir_ratio(const double *points, const void *par,
                            double *weights);

void gen_reference_points_wrapper(const void *par, double *points);

void calc_dummy_ratio(const double * /* unused */, 
                      const void * /* unused */, 
                      double *weights);


char *print_quantiles(const double *quantiles,
                      size_t n_quantiles,
                      const char *label_string,
                      struct locus_sampling *ls, 
                      char *out_buffer);

char *print_discrete_comp(PileupSummary *locus,
                          const char *sample_label,
                          double *discrete_values,
                          size_t n_discrete_values,
                          size_t n_discrete_points_to_print,
                          double min_value_to_print,
                          char *out_buf);

void init_sample_attributes(const char *jpd_file,
                            const char *sample_label,
                            const char *pileup_file,
                            struct sample_attributes *s);

/* initialize the locus in 'ls' defined by sd->current */
void init_pileup_locus(const struct nucleotide_stats *stats,
                       size_t min_quality_score,
                       locus_sampling *ls);


void refresh_locus(const struct nucleotide_stats *stats,
                   size_t min_quality_score,
                   locus_sampling *ls);

thread_queue_worker_t comp_worker;

#endif // _COMP_WORKER_H
