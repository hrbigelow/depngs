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

struct eval_dist_matrix
{
    double * points;
    size_t num_points;
    size_t * dist_index;
    size_t num_distances;
    double * distances;
    double inclusion_threshold; // points with a probability less than this will be ignored
    size_t max_points_to_print;
    double min_value_to_print;

    eval_dist_matrix();
};


enum eval_strategy_t
    {
        LATTICE_ONLY,
        SAMPLING_ONLY,
        LATTICE_THEN_SAMPLING
    };

struct dist_worker_input
{
    posterior_wrapper ** worker; // generally encapsulates everything needed to compute outputs
    size_t thread_num;
    size_t num_samples;
    size_t num_sample_pairs;
    size_t num_sample_point_pairings; 
    double log2_posterior_threshold; // threshold for including a discrete sample point in the weighted pairs.

    double * dist_quantiles;
    size_t num_dist_quantiles;

    double * comp_quantiles;
    size_t num_comp_quantiles;

    float min_high_conf_dist; // minimum mutational distance at high confidence to report change
    gsl_rng * randgen;


    std::vector<char *>::iterator * beg; // ranges for all samples
    std::vector<char *>::iterator * end;
    bool * is_next; // set of flags, one for each sample, defining whether the genomic position of the

    char * out_dist; // points to the next place to write output distances
    char * out_comp; // next place to write output compositions using the sampler
    char * out_discomp; // next place to write output compositions using discrete evaluation
    char * out_vcf; // next place to write output vcf lines

    size_t * pair_sample1; // defines the parsed set of sample pairs to compare
    size_t * pair_sample2;

    std::map<char const*, size_t, ltstr> * contig_order;
    less_locus_position less_locus;
    equal_locus_position equal_locus;

    dist_worker_input(size_t thread_num,
                      size_t num_samples,
                      size_t num_sample_pairs,
                      size_t num_sample_point_pairings,
                      double * dist_quantiles,
                      size_t num_dist_quantiles,
                      double * comp_quantiles,
                      size_t num_comp_quantiles,
                      eval_dist_matrix * lattice,
                      double sampling_fallback_threshold,
                      eval_strategy_t eval_strategy,
                      /* bool use_discrete, */
                      /* bool use_sampling, */
                      char * out_dist,
                      char * out_comp,
                      char * out_discomp,
                      char * out_vcf,
                      size_t * pair_sample1,
                      size_t * pair_sample2,
                      std::map<char const*, size_t, ltstr> * contig_order);

    dist_worker_input();

    ~dist_worker_input();

};


struct sample_details
{
    PileupSummary *locus;
    bool is_next;
    double *sample_points;
    int num_sample_points;
    sampling_method samp_method;
    size_t autocor_offset;
    std::vector<char *>::iterator current;
    char algorithm_used[10];
};


size_t distance_quantiles_locus_bytes(size_t num_quantiles);


void * dist_worker(void * args);


#endif // DIST_WORKER_H
