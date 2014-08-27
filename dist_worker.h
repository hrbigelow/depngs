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
    size_t num_samples;
    size_t num_sample_pairs;
    size_t num_sample_point_pairings; 

    double *dist_quantiles;
    size_t num_dist_quantiles;

    double *comp_quantiles;
    size_t num_comp_quantiles;

    float min_high_conf_dist; // minimum mutational distance at high confidence to report change
    size_t prelim_num_points;
    double prelim_quantile;
    size_t final_num_points;

    gsl_rng *randgen;
    int print_pileup_fields; // if 0, do not print extra pileup fields

    std::vector<char *>::iterator *beg; // ranges for all samples
    std::vector<char *>::iterator *end;
    bool *is_next; // set of flags, one for each sample, defining whether the genomic position of the

    char *out_dist; // points to the next place to write output distances
    char *out_comp; // next place to write output compositions using the sampler
    char *out_vcf; // next place to write output vcf lines
    char *out_indel_dist;

    size_t *pair_sample1, *pair_sample2; // defines the parsed set of sample pairs to compare

    std::map<char const*, size_t, ltstr> *contig_order;
    less_locus_position less_locus;
    equal_locus_position equal_locus;

    dist_worker_input(size_t thread_num,
                      size_t num_samples,
                      size_t num_sample_pairs,
                      size_t num_sample_point_pairings,
                      double *dist_quantiles,
                      size_t num_dist_quantiles,
                      double *comp_quantiles,
                      size_t num_comp_quantiles,
                      float min_high_conf_dist,
                      size_t prelim_num_points,
                      float prelim_quantile,
                      size_t final_num_points,
                      int print_pileup_fields,
                      char *out_dist,
                      char *out_comp,
                      char *out_vcf,
                      char *out_indel_dist,
                      size_t *pair_sample1,
                      size_t *pair_sample2,
                      std::map<char const*, size_t, ltstr> *contig_order);

    dist_worker_input();

    ~dist_worker_input();

};

size_t distance_quantiles_locus_bytes(size_t num_quantiles);


void *dist_worker(void *args);


#endif // DIST_WORKER_H
