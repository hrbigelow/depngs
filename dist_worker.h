#ifndef DIST_WORKER_H
#define DIST_WORKER_H

#include <map>
#include <vector>
#include <cstdlib>
#include <cstring>

struct posterior_wrapper;

struct ltstr
{
    bool operator()(const char* s1, const char* s2) const
    {
        return strcmp(s1, s2) < 0;
    }
};

struct less_locus_position
{
    std::map<char const*, size_t, ltstr> * contig_order;
    bool operator()(char * locus_line1, char * locus_line2);
};

struct equal_locus_position
{
    std::map<char const*, size_t, ltstr> * contig_order;
    bool operator()(char * locus_line1, char * locus_line2);
};


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
};


struct dist_worker_input
{
    posterior_wrapper ** worker; // generally encapsulates everything needed to compute outputs
    size_t num_samples;
    size_t num_sample_pairs;
    size_t num_sample_point_pairings; 
    double log2_posterior_threshold; // threshold for including a discrete sample point in the weighted pairs.

    double * dist_quantiles;
    size_t num_dist_quantiles;

    double * comp_quantiles;
    size_t num_comp_quantiles;

    eval_dist_matrix * lattice;

    double sampling_fallback_threshold;

    bool use_discrete;
    bool use_sampling;

    std::vector<char *>::iterator * beg; // ranges for all samples
    std::vector<char *>::iterator * end;
    bool * is_next; // set of flags, one for each sample, defining whether the genomic position of the

    char * out_dist; // points to the next place to write output distances
    char * out_comp; // next place to write output compositions using the sampler
    char * out_discomp; // next place to write output compositions using discrete evaluation

    size_t * pair_sample1; // defines the parsed set of sample pairs to compare
    size_t * pair_sample2;

    std::map<char const*, size_t, ltstr> * contig_order;
    less_locus_position less_locus;
    equal_locus_position equal_locus;

    dist_worker_input(size_t num_samples,
                      size_t num_sample_pairs,
                      size_t num_sample_point_pairings,
                      double * dist_quantiles,
                      size_t num_dist_quantiles,
                      double * comp_quantiles,
                      size_t num_comp_quantiles,
                      eval_dist_matrix * lattice,
                      double sampling_fallback_threshold,
                      bool use_discrete,
                      bool use_sampling,
                      char * out_dist,
                      char * out_comp,
                      char * out_discomp,
                      size_t * pair_sample1,
                      size_t * pair_sample2,
                      std::map<char const*, size_t, ltstr> * contig_order);

    dist_worker_input();

    ~dist_worker_input();

};

size_t distance_quantiles_locus_bytes(size_t num_quantiles);


void * dist_worker(void * args);


#endif // DIST_WORKER_H
