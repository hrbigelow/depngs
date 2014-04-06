#include "anomaly_tools.h"

#include <cassert>
#include <numeric>

#include "nucleotide_stats.h"
#include "integrands.h"
#include "error_estimate.h"
#include "stats_tools.h"

#include <gsl/gsl_math.h>


 // given P(basecall, quality, strand | founder_base),
 // first compute P(basecall, quality | founder_base, strand = S)
 // then 
void
normalize_strand_marginal(size_t strand_to_marginalize,
                          double const* p_bqs_given_f,
                          double * p_bqs_new,
                          size_t * index_map,
                          size_t num_data)
{
    
    // iterate over the initial distribution, generating a
    // distribution normalized by strand.
    size_t S = strand_to_marginalize;
    size_t notS = 1 - S;

    char basecall;
    size_t quality;
    size_t strand;

    double strand_marg[2] = {0.0, 0.0};
    for (size_t i = 0; i != num_data; ++i)
    {
        Nucleotide::decode(index_map[i], &basecall, &quality, &strand);
        strand_marg[strand] += p_bqs_given_f[i];
    }

    // in-place normalization
    normalize(strand_marg, 2, strand_marg);

    double flow[2];
    flow[S] = 1.0 / strand_marg[S];
    flow[notS] = 0.0;

    assert(! isnan(flow[S]));

    for (size_t i = 0; i != num_data; ++i)
    {
        Nucleotide::decode(index_map[i], &basecall, &quality, &strand);
        p_bqs_new[i] = p_bqs_given_f[i] * flow[strand];
    }


}


/*
  computes D. A. Williams' adjustment factor for improving the likelihood ratio
  see Wikipedia's entry on the "Multinomial Test"
  Williams, D. A. (1976). "Improved Likelihood Ratio Tests for Complete Contingency Tables".
  Biometrika 63: 33-37. http://dx.doi.org/10.1093/biomet/63.1.33
  
  assumes any pi factors equalling zero are not involved, so doesn't compute their reciprocal
*/
double williams_moment_match_ratio(NucleotideStats const& stats,
                                   size_t N,
                                   double const* evaluation_point)
{
    size_t k = 0;
    
    double sum_pi = 0.0;
    double sum_pi_reciprocal = 0.0;
    for (size_t i = 0; i != Nucleotide::num_bqs; ++i)
    {
        double const* mp = evaluation_point;
            
        double pi_i =
            stats.founder_base_likelihood[0][i] * mp[0]
            + stats.founder_base_likelihood[1][i] * mp[1]
            + stats.founder_base_likelihood[2][i] * mp[2]
            + stats.founder_base_likelihood[3][i] * mp[3];
            
        sum_pi += pi_i;

        if (pi_i != 0)
        {
            ++k;
            sum_pi_reciprocal += 1.0 / pi_i;
        }
    }

    assert(gsl_fcmp(sum_pi, 1.0, 1e-10) == 0);

    double ratio = 1.0 + (sum_pi_reciprocal - 1.0) / (6.0 * N * (k - 1));
        
    return ratio;
}


/*
  The idea here is to separately collect the positive and negative
  stranded locus reads data, and run it through the model.

  The model must be parameterized in a strand-specific way, that
  is, it must be blind to strand.

  Finally, in order to get an anomaly score, we run the data
  through two models: The model parameterized by global data (data
  from all loci), and secondly, the model parameterized by just
  the locus data for the locus in question.

  Because the anomaly score attempts to find sequence-context
  specific effects, it makes sense to run it separately on
  positive and negative stranded data.

*/
/*
  double strand_locus_anomaly_score(Posterior & posterior,
  JPD_DATA const& global_counts,
  LocusSummary const& full_locus,
  NucleotideReader const* data_reader,
  char strand,
  bool verbose)
  {

  double const mode_tolerance = 1e-10;
  double const max_modefinding_iterations = 3000;

  LocusSummary strand_locus =
  data_reader->locus_data_by_strand(full_locus, strand);

                                         
  double strand_marginal[] = { 0, 0 };
  if (strand == '+')
  {
  strand_marginal[0] = 1;
  }
  else
  {
  strand_marginal[1] = 1;
  }

  JPD_DATA strand_global_counts = 
  data_reader->normalize_strand_marginal(strand_marginal, global_counts);

  NucleotideStats & post_params = * posterior.ee->get_model_params();
  //post_params.initialize(global_counts);
  post_params.initialize(strand_global_counts);

  JPD_DATA strand_locus_counts = 
  post_params.make_per_locus_stats(strand_locus);
    
  //posterior.ee->set_locus_data(& full_locus);
  posterior.ee->set_locus_data(& strand_locus);

  double initial_point[] = { 0.25, 0.25, 0.25, 0.25 };
  posterior.initialize(mode_tolerance, max_modefinding_iterations, initial_point, verbose);

  // double williams_strand_ratio = 
  //     williams_moment_match_ratio(post_params, strand_locus.read_depth, 
  //                                 posterior.mode_point);


  double strand_global_log_peak = posterior.log_pdf(posterior.mode_point);

  post_params.initialize(strand_locus_counts);

  posterior.initialize(mode_tolerance, max_modefinding_iterations, initial_point, verbose);
  double strand_locus_log_peak = posterior.log_pdf(posterior.mode_point);


  double strand_anomaly_score = 
  -2.0 * (strand_global_log_peak - strand_locus_log_peak)
  / strand_locus.read_depth;
    
  // / williams_strand_ratio;

  return strand_anomaly_score;

  }
*/


// evaluate the posterior peak value two ways:
// 1. parameterize the model with locus-derived counts
// 2. parameterize the model with global counts
// Both model parameters and locus data are the subset pertaining to a
// single chosen strand.  This is to eliminate any anomaly due to stranded behavior.
/*
double strand_locus_anomaly_score(Posterior & posterior,
                                  double const* cpd_fbqs_global,
                                  packed_counts * locus_counts,
                                  size_t which_strand,
                                  bool verbose)
{

    double const mode_tolerance = 1e-10;
    double const max_modefinding_iterations = 3000;

    LocusSummary strand_locus =
        data_reader->locus_data_by_strand(full_locus, strand);

    // global_counts is the Phred-inferred C(basecall, quality, strand, founder_base) distribution
    // normalize_strand_marginal(strand, global_stats.jpd_buffer,
    //                           norm_jpd_buffer,
                              

    JPD_DATA strand_global_counts = 
        data_reader->normalize_strand_marginal(strand_marginal, global_counts);

    ErrorEstimate * model = posterior.ee;
    NucleotideStats * post_params = model->model_params;

    post_params->initialize(strand_global_counts);
    JPD_DATA strand_locus_counts = post_params->make_per_locus_stats(strand_locus);
    model->locus_data = strand_counts;

    double initial_point[] = { 0.25, 0.25, 0.25, 0.25 };
    posterior.initialize(mode_tolerance, max_modefinding_iterations, initial_point, verbose);

    double strand_global_log_peak = posterior.log_pdf(posterior.mode_point);

    post_params->initialize(strand_locus_counts);

    posterior.initialize(mode_tolerance, max_modefinding_iterations, initial_point, verbose);
    double strand_locus_log_peak = posterior.log_pdf(posterior.mode_point);


    double strand_anomaly_score = 
        -2.0 * (strand_global_log_peak - strand_locus_log_peak)
        / strand_locus.read_depth;
    
    // / williams_strand_ratio;

    return strand_anomaly_score;
}
*/


/*
//compute relative entropy between P_global(b,q,s) and P_local(b,q,s)
//P_global(b,q,s) is computed from P_global(b,q,s,f) using a founder base (f)
//composition given as argument (typically the median)
double relative_entropy_anomaly(JPD_DATA const& global_counts,
LocusSummary const& locus_summary,
double est_composition[4])
{
    
NucleotideStats params(global_counts.size());
params.initialize(global_counts);

double * p_global = new double[params.num_distinct_data];
double * p_local = new double[params.num_distinct_data];

size_t const D = params.num_distinct_data;

//compute P_global(b,q,s)
for (size_t d = 0; d != D; ++d)
{
p_global[d] = 
params.founder_base_likelihood[0][d] * est_composition[0]
+ params.founder_base_likelihood[1][d] * est_composition[1]
+ params.founder_base_likelihood[2][d] * est_composition[2]
+ params.founder_base_likelihood[3][d] * est_composition[3];
}

normalize(p_global, D, p_global);

//compute P_local(b,q,s), simply marginalize P_local(b,q,s,f)
JPD_DATA local_counts = params.make_per_locus_stats(locus_summary);

params.initialize(local_counts);
for (size_t d = 0; d != D; ++d)
{
p_local[d] =
params.complete_jpd[0][d]
+ params.complete_jpd[1][d]
+ params.complete_jpd[2][d]
+ params.complete_jpd[3][d];
}
normalize(p_local, D, p_local);

//    for (size_t d = 0; d != D; ++d)
//    {
//        printf("%s\t%Zu\t%f\t%f\n", params.index_mapping[d].c_str(),
//               d, p_local[d], p_global[d]);
//    }

double kl_divergence = kullback_leibler_divergence(p_local, p_global, D);

delete p_global;
delete p_local;

return kl_divergence;
}
*/

 //compute relative entropy between P_global(b,q,s) and P_local(b,q,s)
 //P_global(b,q,s) is computed from P_global(b,q,s,f) using a founder base (f)
 //composition given as argument (typically the median)
double relative_entropy_anomaly(double const* cpd_fbqs_global,
                                packed_counts * locus_stats,
                                double est_composition[4])
{
    size_t const D = locus_stats->num_data;

    double * p_global = new double[D];
    double * p_local = new double[D];

    //compute P_global(b,q,s)
    double const* g = cpd_fbqs_global;
    double const* g_end = g + Nucleotide::num_bqs * 4;
    size_t d = 0;
    for (; g != g_end; g += 4, ++d)
    {
        p_global[d] = std::inner_product(g, g + 4, est_composition, 0);
        // params.founder_base_likelihood[0][d] * est_composition[0]
        // + params.founder_base_likelihood[1][d] * est_composition[1]
        // + params.founder_base_likelihood[2][d] * est_composition[2]
        // + params.founder_base_likelihood[3][d] * est_composition[3];
    }

    normalize(p_global, D, p_global);

    //compute P_local(b,q,s), simply marginalize P_local(b,q,s,f)

    double * ls = locus_stats->fbqs_cpd;
    double * ls_end = locus_stats->fbqs_cpd + locus_stats->num_data;
    for (; ls != ls_end; ls += 4)
    {
        p_local[d] = std::accumulate(ls, ls + 4, 0);
    }
    normalize(p_local, D, p_local);

    double kl_divergence = kullback_leibler_divergence(p_local, p_global, D);

    delete p_global;
    delete p_local;

    return kl_divergence;
}

/*
double locus_anomaly_score(Posterior & posterior,
                           JPD_DATA const& global_counts,
                           LocusSummary const& full_locus,
                           NucleotideReader const* data_reader,
                           bool verbose)
{

    double const mode_tolerance = 1e-10;
    double const max_modefinding_iterations = 3000;

    NucleotideStats & post_params = * posterior.ee->get_model_params();
    post_params.initialize(global_counts);

    JPD_DATA locus_counts = post_params.make_per_locus_stats(full_locus);
    
    posterior.ee->set_locus_data(& full_locus);

    double initial_point[] = { 0.25, 0.25, 0.25, 0.25 };
    posterior.initialize(mode_tolerance, max_modefinding_iterations, 
                         initial_point, verbose);

    // double williams_strand_ratio = 
    //     williams_moment_match_ratio(post_params, strand_locus.read_depth, 
    //                                 posterior.mode_point);

    double global_log_peak = posterior.log_pdf(posterior.mode_point);

    post_params.initialize(locus_counts);

    posterior.initialize(mode_tolerance, max_modefinding_iterations, 
                         initial_point, verbose);

    double locus_log_peak = posterior.log_pdf(posterior.mode_point);

    double anomaly_score = 
        -2.0 * (global_log_peak - locus_log_peak)
        / full_locus.read_depth;
    
    // / williams_strand_ratio;

    return anomaly_score;

}
*/


/*
  ErrorEstimate simulation_model;
  simulation_model.set_model_params(& global_params);
  simulation_model.set_composition_prior_alphas(prior_alphas);

  double * log_simulated_modes = new double[nsim_loci];
  for (size_t si = 0; si != nsim_loci; ++si)
  {
  LocusSummary simulated_locus = 
  sample_locus_from_stats(rand_gen, global_params, 
  posterior.mode_point,
  locus.read_depth);

  simulation_model.set_locus_data(&simulated_locus);
  Posterior simulation_posterior(&simulation_model, may_underflow, 
  full_ndim);

  //              simulation_posterior.initialize(mode_tolerance, 
  //                                              max_modefinding_iterations, 0);

  log_simulated_modes[si] = 
  simulation_posterior.log_pdf(posterior.mode_point);

  }
  //calculate mean and sd of log values and regular values
  double mean_log_simulated = 
  gsl_stats_mean(log_simulated_modes, 1, nsim_loci);

  double sd_log_simulated =
  gsl_stats_sd_m(log_simulated_modes, 1, nsim_loci, 
  mean_log_simulated);

  double anomaly_score = 
  (log_this_mode - mean_log_simulated)
  / sd_log_simulated;

  delete log_simulated_modes;



  //set local model posterior
  //         NucleotideStats local_params = global_params.make_per_locus_stats(locus);
  //         ErrorEstimate local_model;
  //         local_model.set_locus_data(& locus);
  //         local_model.set_model_params(& local_params);
  //         Posterior this_locus_posterior(&local_model, may_underflow, truncated_ndim);
  //         this_locus_posterior.initialize(mode_tolerance, max_modefinding_iterations, 0);


  //         double locus_log_peak = this_locus_posterior.log_pdf(this_locus_posterior.mode_point);
  //         double main_log_peak = posterior.log_pdf(posterior.mode_point);

  //assert(locus_log_peak >= main_log_peak);

  //         double anomaly_score = main_log_peak - locus_log_peak;
  */
