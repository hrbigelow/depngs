#ifndef _ERROR_ESTIMATE_H
#define _ERROR_ESTIMATE_H

#include <cstddef>
#include <map>

struct packed_counts;

/*
Estimate the distribution of base compositions for a given sequenced sample.

See full_posterior_model.doc

P(m | (b1,q1),...,(bn,qn)) as
prod { P((bi,qi)|m) * P(m) / P(bi,qi) } as
prod { P((bi,qi)|m) * P(m) } / prod { P(bi,qi) }

(in log form)

sum { log(P((bi,qi)|m) * P(m)) } - sum { log(P(bi,qi)) }
*/
                                  //enum DNAStrand { POS_STRAND, NEG_STRAND };

class NucleotideStats;
class LocusSummary;


class ErrorEstimate {

 public:
    static char const* nucleotides;
    static int const base_to_index[];

 private:

    //for expanding to the unit hypercube
    double expansion_rows[3][3];

    //for contracting from the unit hypercube
    double contraction_rows[3][3];


 public:

    ErrorEstimate();
    ~ErrorEstimate();

    packed_counts * locus_data;
    NucleotideStats * model_params;

    double composition_prior_alphas[4];
    
    void set_composition_prior_alphas(double const* alphas);

    void set_discrete_prior_dist(double const* prior_points_flat,
                                 double const* prior_dist, 
                                 size_t num_priors);

    void log_likelihood_gradient(double const* sample_composition,
                                 double * gradient) const;
    
    double log_discrete_prior(size_t sample_point_index) const;

    double log_dirichlet_prior(double const* sample_composition) const;

    double log_likelihood(double const* sample_composition) const;

        
    double ScaledPosterior(double const* sample_composition,
                           double log_scaling_factor) const;
        

    void expand_to_hypercube(double const x[3], double * expanded) const;
    void contract_from_hypercube(double const x[3], double * condensed) const;

    size_t find_mode_point(double min_step_size,
                           size_t max_iterations,
                           double const* initial_point,
                           bool * on_zero_boundary,
                           bool verbose,
                           double * mode_point) const;

        
};

double log2_likelihood(ErrorEstimate * ee, double const* comp);

//tells whether x is within the unit hypercube
bool within_hypercube(double const x[3]);

//tells whether x is a coordinate within the corner pyramid
bool within_pyramid(double const x[3]);

#endif // _ERROR_ESTIMATE_H
