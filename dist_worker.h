#ifndef DIST_WORKER_H
#define DIST_WORKER_H

#include <map>
#include <vector>
#include <cstdlib>
#include <cstring>

#include <gsl/gsl_randist.h>

#include "locus_comp.h"

struct posterior_wrapper;

class PileupSummary;


struct dist_worker_input
{
    posterior_wrapper **worker; // generally encapsulates everything needed to compute outputs
    size_t thread_num;
    size_t n_samples;
    size_t n_sample_pairs;
    size_t n_sample_point_pairings; 

    double *dist_quantiles;
    size_t n_dist_quantiles;

    double *comp_quantiles;
    size_t n_comp_quantiles;

    float min_high_conf_dist; // minimum mutational distance at high confidence to report change
    size_t prelim_n_points;
    double prelim_quantile;
    size_t final_n_points;

    gsl_rng *randgen;
    int print_pileup_fields; // if 0, do not print extra pileup fields

    /* std::vector<char *>::iterator *beg */
    /* std::vector<char *>::iterator *end; */
    bool *is_next; // set of flags, one for each sample, defining whether the genomic position of the

    /* indices of the output buffers.  if -1, do not perform this
       function */
    int do_dist, do_comp, do_indel, do_vcf;

    /* defines the parsed set of sample pairs to compare */
    size_t *pair_sample1, *pair_sample2; 

    less_locus_position less_locus;
    equal_locus_position equal_locus;

    dist_worker_input(size_t thread_num,
                      size_t n_samples,
                      size_t n_sample_pairs,
                      size_t n_sample_point_pairings,
                      double *dist_quantiles,
                      size_t n_dist_quantiles,
                      double *comp_quantiles,
                      size_t n_comp_quantiles,
                      float min_high_conf_dist,
                      size_t prelim_n_points,
                      float prelim_quantile,
                      size_t final_n_points,
                      int print_pileup_fields,
                      int do_dist,
                      int do_comp,
                      int do_indel,
                      int do_vcf,
                      size_t *pair_sample1,
                      size_t *pair_sample2);

    dist_worker_input();

    ~dist_worker_input();

};


struct dist_worker_offload_par {
    FILE *dist_fh, *comp_fh, *indel_fh, *vcf_fh;
};

/* conforms to thread_queue_worker_t */
void dist_worker(void *par,
                 const struct managed_buf *in_bufs,
                 struct managed_buf *out_bufs);

/* conforms to thread_queue_offload_t */
void dist_offload(void *par, const struct managed_buf *bufs);

#endif // DIST_WORKER_H
