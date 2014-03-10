#include <cstddef>
#include <numeric>

#include "nucleotide_stats.h"
#include "error_estimate.h"
#include "dirichlet.h"
#include "metropolis.h"
#include "pileup_tools.h"
#include "stats_tools.h"
#include "slice_sampling.h"

#include "comp_functor.h"

/*
  input: the null-terminated pileup line
  side-effects: populate the appropriate line buffer with the output
 */



posterior_wrapper::posterior_wrapper(char const* jpd_data_params_file,
                                     double * prior_alphas,
                                     size_t min_quality_score,
                                     double * quantiles,
                                     size_t num_quantiles,
                                     char const* label_string,
                                     FILE * cdfs_output_fh,
                                     pthread_mutex_t * file_writing_mutex,
                                     bool compute_anomaly) :
    mode_tolerance(1e-60),
    max_modefinding_iterations(3000),
    max_tuning_iterations(10),
    tuning_num_points(1000),
    final_num_points(10000),
    autocor_max_offset(6),
    autocor_max_value(6),
    initial_autocor_offset(30),
    target_autocor_offset(6),
    num_bits_per_dim(62),
    is_log_integrand(true),
    initial_sampling_range(62 * 3),
    may_underflow(true),
    use_independence_chain_mh(true),
    verbose(false),
    min_quality_score(min_quality_score),
    num_quantiles(num_quantiles),
    cdfs_output_fh(cdfs_output_fh),
    file_writing_mutex(file_writing_mutex),
    compute_anomaly(compute_anomaly)
{
    this->label_string = new char[strlen(label_string) + 1];
    strcpy(this->label_string, label_string);

    this->quantiles = new double[this->num_quantiles];
    std::copy(quantiles, quantiles + num_quantiles, this->quantiles);

    this->initial_point = { 0.25, 0.25, 0.25, 0.25 };
    size_t const full_ndim = 4;

    this->sample_points_buf = new double[this->final_num_points * 4];

    this->prior_alpha0 = std::accumulate(prior_alphas, prior_alphas + 4, 0.0);
    
    this->params = new NucleotideStats();
    this->params->initialize(jpd_data_params_file);
    this->model = new ErrorEstimate();
    this->model->set_composition_prior_alphas(prior_alphas);
    this->model->model_params = this->params;
    this->prior = new Dirichlet(full_ndim, this->may_underflow);
    this->posterior = new Posterior(this->model, this->may_underflow, full_ndim);
    this->sampler = new Metropolis(this->posterior, this->prior, full_ndim, 
                                   this->use_independence_chain_mh, this->final_num_points);

    // this->sample_points_sortable.resize(this->final_num_points);
    // for (size_t i = 0; i != this->final_num_points; ++i)
    // {
    //     this->sample_points_sortable[i] = this->sampler->sample_points + (i * full_ndim);
    // }

}



posterior_wrapper::~posterior_wrapper()
{
    delete this->label_string;
    delete this->quantiles;

    delete this->sample_points_buf;
    delete this->sampler;
    delete this->posterior;
    delete this->prior;
    delete this->model;
    delete this->params;
}



void posterior_wrapper::process_line_comp(char const* pileup_line,
                                          char * out_buffer)
{

    PileupSummary locus(0);
    locus.load_line(pileup_line);
    locus.parse(this->min_quality_score);
    this->params->pack(& locus.counts);

    //divide locus data to plus and minus-strand data
    // double pos_anomaly_score =
    //     strand_locus_anomaly_score(posterior, params.,
    //                                locus, data_reader, '+', verbose);
            
    // double neg_anomaly_score =
    //     strand_locus_anomaly_score(posterior, global_counts,
    //                                locus, data_reader, '-', verbose);
            
            
    // double full_anomaly_score =
    //     locus_anomaly_score(posterior, global_counts,
    //                         locus, data_reader, verbose);


    this->posterior->model()->locus_data = & locus.counts;
    this->posterior->initialize(this->mode_tolerance, this->max_modefinding_iterations, 
                                initial_point, verbose);
    this->sampler->set_current_point(this->posterior->mode_point);
    this->prior->set_alpha0(locus.read_depth + this->prior_alpha0);

    if (this->prior->get_alpha0() > 1)
    {
        this->prior->set_alphas_from_mode_or_bound
            (this->posterior->mode_point,
             this->posterior->ee->composition_prior_alphas,
             this->posterior->zero_boundary);
    }
    else
    {
        this->prior->update(this->posterior->ee->composition_prior_alphas);
    }

    //metropolis hastings
    double proposal_mean, proposal_variance;
    size_t cumul_autocor_offset;
    size_t best_autocor_offset = this->initial_autocor_offset;
    size_t const full_ndim = 4;
    size_t const truncated_ndim = 3;
    double estimated_mean[full_ndim];
    char line_label[1000];
    size_t effective_depth = locus.read_depth;

    char const* dimension_labels[] = { "A", "C", "G", "T" };

    for (size_t iter = 0; iter != this->max_tuning_iterations; ++iter)
    {
        cumul_autocor_offset = 1;
        size_t current_autocor_offset = 1;
        for (size_t i = 0; i != 3; ++i)
        {
            //sample more and more thinly, starting from every 1'th
            this->sampler->sample(this->tuning_num_points, 0, cumul_autocor_offset,
                                  &proposal_mean, &proposal_variance);
            
            current_autocor_offset =
                best_autocorrelation_offset(this->sampler->sample_points,
                                            full_ndim, this->tuning_num_points,
                                            this->autocor_max_offset, this->autocor_max_value);

            if (current_autocor_offset == 1)
            {
                break;
            }

            cumul_autocor_offset *= current_autocor_offset;
            if (verbose)
            {
                fprintf(stdout, "MH: current: %Zu, cumulative: %Zu, position: %i\n", 
                        current_autocor_offset, cumul_autocor_offset, 
                        locus.position);
                fflush(stdout);
            }
        }

        //here it doesn't make sense to take all of the samples if
        //the best_autoocr offset isn't 1
        multivariate_mean(this->sampler->sample_points, full_ndim,
                          this->tuning_num_points, estimated_mean);

        if (verbose)
        {
            fprintf(stdout, "MH: %Zu, position: %i, proposal mean: %g, "
                    "proposal variance: %g", cumul_autocor_offset,
                    locus.position, proposal_mean, proposal_variance);
            fprintf(stdout, ", estimated mean:");
            for (size_t d = 0; d != full_ndim; ++d)
            {
                fprintf(stdout, "\t%g", estimated_mean[d]);
            }
            fprintf(stdout, "\n");
            fflush(stdout);
        }

        if (cumul_autocor_offset <= target_autocor_offset)
        {
            break;
        }

        this->prior->set_alphas_from_mean_or_bound(estimated_mean,
                                                   this->posterior->ee->composition_prior_alphas);
        
    }
    
    if (cumul_autocor_offset <= this->target_autocor_offset)
    {
        //metropolis hastings succeeded
        this->sampler->sample(this->final_num_points, 0, cumul_autocor_offset,
                             &proposal_mean, &proposal_variance);

        sprintf(line_label, "%s\t%s\t%s\t%i\t%c\t%Zu\t%Zu"
                // "\t%5.5lf\t%5.5lf\t%5.5lf"
                ,
                label_string, 
                "MH", locus.reference, 
                locus.position, 
                locus.reference_base, 
                locus.read_depth, 
                effective_depth
                // ,
                // full_anomaly_score,
                // pos_anomaly_score,
                // neg_anomaly_score
                );

        //             //approximate marginal modes (expensive!)
        //             std::fill(marginal_modes, marginal_modes + 4, -1.0);
        //             if (marginal_mode_num_points > 0)
        //             {
        //                 for (size_t d = 0; d != full_ndim; ++d)
        //                 {
        //                     marginal_modes[d] = 
        //                         window_averaged_mode(& sample_points_sortable, d, 
        //                                              marginal_mode_num_points);
        //                 }
        //             }

        print_quantiles(out_buffer,
                        this->sampler->sample_points,
                        this->final_num_points,
                        this->posterior->mode_point,
                        line_label, dimension_labels, "+", quantiles,
                        num_quantiles, full_ndim);

        if (cdfs_output_fh != NULL)
        {
            pthread_mutex_lock(file_writing_mutex);
            print_numerical_cdfs(cdfs_output_fh, 
                                 line_label, 
                                 this->sampler->sample_points,
                                 this->final_num_points, 
                                 full_ndim);

            pthread_mutex_unlock(file_writing_mutex);
        }
        
        return;
    }

    if (verbose)
    {
        fprintf(stderr, 
                "%s %i (depth %Zu) locus skipped"
                ", Metropolis Hastings failed.\n",
                locus.reference, locus.position, effective_depth);
    }

    //metropolis hastings failed.  try slice sampling
    SliceSampling slice_sampler(truncated_ndim, num_bits_per_dim, is_log_integrand, 1);

    double * tuning_sample_points = new double[tuning_num_points * truncated_ndim];

    slice_sampler.Initialize();

    slice_sampler.sample(this->posterior, 
                         this->posterior->mode_point,
                         initial_sampling_range, 1, tuning_sample_points,
                         tuning_num_points);
        
    best_autocor_offset =
        best_autocorrelation_offset(tuning_sample_points, truncated_ndim, 
                                    tuning_num_points, autocor_max_offset, 
                                    autocor_max_value);
    if (verbose)
    {
        fprintf(stdout, "SS: %Zu, position: %i\n", 
                best_autocor_offset, locus.position);
    }
        
    delete tuning_sample_points;


    if (best_autocor_offset <= target_autocor_offset)
    {
        //slice sampling succeeded
        double * slice_sample_points = new double[this->final_num_points * full_ndim];

        slice_sampler.sample(this->posterior, 
                             this->posterior->mode_point,
                             initial_sampling_range, 
                             best_autocor_offset,
                             slice_sample_points,
                             this->final_num_points);

        sprintf(line_label, "%s\t%s\t%s\t%i\t%c\t%Zu\t%Zu"
                // "\t%5.5lf\t%5.5lf\t%5.5lf"
                ,
                label_string, 
                "SS", locus.reference, 
                locus.position, 
                locus.reference_base, 
                locus.read_depth, 
                effective_depth
                // ,
                // full_anomaly_score,
                // pos_anomaly_score, 
                // neg_anomaly_score
                );
            
        //approximate marginal modes (expensive!)
        //             std::fill(marginal_modes, marginal_modes + 4, -1.0);
        //             if (marginal_mode_num_points > 0)
        //             {
        //                 for (size_t d = 0; d != full_ndim; ++d)
        //                 {
        //                     marginal_modes[d] = 
        //                         window_averaged_mode(& slice_sample_points_sortable, d, 
        //                                              marginal_mode_num_points);
        //                 }
        //             }
            
        print_quantiles(out_buffer, 
                        slice_sample_points, 
                        this->final_num_points,
                        this->posterior->mode_point, 
                        line_label, 
                        dimension_labels, "+", quantiles, 
                        num_quantiles,
                        truncated_ndim + 1);


        if (cdfs_output_fh != NULL)
        {
            pthread_mutex_lock(file_writing_mutex);
            print_numerical_cdfs(cdfs_output_fh, line_label, slice_sample_points, this->final_num_points, full_ndim);
            pthread_mutex_unlock(file_writing_mutex);
        }

        delete slice_sample_points;
    }
}



void posterior_wrapper::process_line_mode(char const* pileup_line,
                                          char * out_buffer)
{

    PileupSummary locus(0);
    locus.load_line(pileup_line);
    locus.parse(this->min_quality_score);
    this->params->pack(& locus.counts);

    //divide locus data to plus and minus-strand data
    // double pos_anomaly_score =
    //     strand_locus_anomaly_score(posterior, params.,
    //                                locus, data_reader, '+', verbose);
            
    // double neg_anomaly_score =
    //     strand_locus_anomaly_score(posterior, global_counts,
    //                                locus, data_reader, '-', verbose);
            
            
    // double full_anomaly_score =
    //     locus_anomaly_score(posterior, global_counts,
    //                         locus, data_reader, verbose);


    size_t effective_depth = locus.read_depth;
    this->posterior->model()->locus_data = & locus.counts;
    this->posterior->initialize(this->mode_tolerance, this->max_modefinding_iterations, 
                                initial_point, verbose);

    out_buffer += 
        sprintf(out_buffer, 
                "%s\t%s\t%i\t%c\t%Zu\t%Zu\t%5.5f\t%5.5f\t%5.5f\t%5.5f",
                label_string, 
                locus.reference, 
                locus.position, 
                locus.reference_base, 
                locus.read_depth,
                effective_depth,
                this->posterior->mode_point[0],
                this->posterior->mode_point[1],
                this->posterior->mode_point[2],
                this->posterior->mode_point[3]
                );
    
    if (compute_anomaly)
    {
        // double pos_anomaly_score =
        //     strand_locus_anomaly_score(posterior, global_counts,
        //                                locus, data_reader, '+', verbose);
        
        // double neg_anomaly_score =
        //     strand_locus_anomaly_score(posterior, global_counts,
        //                                locus, data_reader, '-', verbose);
        // double full_anomaly_score =
        //     relative_entropy_anomaly(params.cpd_buffer, & locus.counts, posterior.mode_point);
        
        // sprintf(out_buffer,
        //         "\t%5.5lf\t%5.5lf\t%5.5lf\n",
        //         full_anomaly_score,
        //         pos_anomaly_score,
        //         neg_anomaly_score);
    }
    else
    {
        sprintf(out_buffer, "\n");
    }
    
}


void * comp_worker(void * args)
{
    wrapper_input * input = static_cast<wrapper_input *>(args);
    std::vector<char *>::iterator it;
    char ** out = input->out_start;
    for (it = input->beg; it != input->end; ++it)
    {
        input->worker->process_line_comp(*it, *out);
        ++out;
    }
    pthread_exit((void*) 0);
}



void * mode_worker(void * args)
{
    wrapper_input * input = static_cast<wrapper_input *>(args);
    std::vector<char *>::iterator it;
    char ** out = input->out_start;
    for (it = input->beg; it != input->end; ++it)
    {
        input->worker->process_line_mode(*it, *out);
        ++out;
    }
    pthread_exit((void*) 0);
}
