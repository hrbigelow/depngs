#include "comp_functor.h"

#include <cstring>
#include <numeric>
#include <algorithm>

#include <gsl/gsl_sf_exp.h>

#include "nucleotide_stats.h"
#include "error_estimate.h"
#include "dirichlet.h"
#include "metropolis.h"
#include "pileup_tools.h"
#include "stats_tools.h"
#include "slice_sampling.h"
#include "error_estimate.h"

// input: the null-terminated pileup line
// side-effects: populate the appropriate line buffer with the output


sample_details::sample_details(void) :
    locus(NULL), is_next(false), sample_points(NULL),
    num_sample_points(0), samp_method(FAILED),
    mode_computed(false),
    autocor_offset(1000),
    current(NULL)
{}


posterior_wrapper::posterior_wrapper(const char *jpd_data_params_file,
                                     double *prior_alphas,
                                     size_t min_quality_score,
                                     double *quantiles,
                                     size_t num_quantiles,
                                     const char *label_string,
                                     FILE *cdfs_output_fh,
                                     pthread_mutex_t *file_writing_mutex,
                                     double gradient_tolerance,
                                     size_t tuning_num_points,
                                     size_t final_num_points,
                                     bool verbose) :
    gradient_tolerance(gradient_tolerance),
    max_modefinding_iterations(3000),
    max_tuning_iterations(10),
    tuning_num_points(tuning_num_points),
    final_num_points(final_num_points),
    autocor_max_offset(6),
    autocor_max_value(6),
    initial_autocor_offset(30),
    target_autocor_offset(6),
    num_bits_per_dim(62),
    is_log_integrand(true),
    initial_sampling_range(62 * 3),
    may_underflow(true),
    use_independence_chain_mh(true),
    verbose(verbose),
    min_quality_score(min_quality_score),
    num_quantiles(num_quantiles),
    cdfs_output_fh(cdfs_output_fh),
    file_writing_mutex(file_writing_mutex)
{
    this->label_string = new char[strlen(label_string) + 1];
    strcpy(this->label_string, label_string);

    this->quantiles = new double[this->num_quantiles];
    std::copy(quantiles, quantiles + num_quantiles, this->quantiles);

    this->initial_point = { 0.25, 0.25, 0.25, 0.25 };
    const size_t full_ndim = 4;
    const size_t truncated_ndim = 3;

    this->prior_alpha0 = std::accumulate(prior_alphas, prior_alphas + 4, 0.0);
    
    this->params = new NucleotideStats();
    this->params->initialize(jpd_data_params_file);
    this->model = new ErrorEstimate();
    this->model->set_composition_prior_alphas(prior_alphas);
    this->model->model_params = this->params;
    this->prior = new Dirichlet(full_ndim, this->may_underflow);
    // this->posterior = new Posterior(this->model, this->may_underflow, full_ndim);
    this->sampler = new Metropolis(this->model, this->prior, full_ndim, 
                                   this->use_independence_chain_mh, this->final_num_points);
    this->slice_sampler = new SliceSampling(truncated_ndim, this->num_bits_per_dim, 1);
}



posterior_wrapper::~posterior_wrapper()
{
    delete this->label_string;
    delete this->quantiles;

    delete this->sampler;
    delete this->slice_sampler;
    // delete this->posterior;
    delete this->prior;
    delete this->model;
    delete this->params;
}



void posterior_wrapper::find_mode(void)
{
    if (this->model->locus_data->num_data == 0)
    {
        memcpy(this->mode_point, NULL_MODE, sizeof(double) * 4);
        this->on_zero_boundary[0] = false;
        this->on_zero_boundary[1] = false;
        this->on_zero_boundary[2] = false;
        this->on_zero_boundary[3] = false;
    }
    else
        this->model->find_mode_point(this->gradient_tolerance,
                                     this->max_modefinding_iterations,
                                     this->initial_point,
                                     this->on_zero_boundary,
                                     this->verbose,
                                     this->mode_point);
}


// tune the posterior, return the cumulative autocorrelation offset
// sample_points_buf is used as a space for writing sample points
// during the tuning process.
size_t posterior_wrapper::tune_mh(PileupSummary *locus, double *sample_points_buf)
{

    this->sampler->set_current_point(this->mode_point);
    this->prior->set_alpha0(locus->read_depth + this->prior_alpha0);

    if (this->prior->get_alpha0() > 1)
    {
        this->prior->set_alphas_from_mode_or_bound
            (this->mode_point,
             this->model->composition_prior_alphas,
             this->on_zero_boundary);
    }
    else
    {
        this->prior->update(this->model->composition_prior_alphas);
    }

    //metropolis hastings
    double proposal_mean, proposal_variance;
    size_t cumul_autocor_offset = 0;
    const size_t full_ndim = 4;
    double estimated_mean[full_ndim];

    
    for (size_t iter = 0; iter != this->max_tuning_iterations; ++iter)
    {
        cumul_autocor_offset = 1;
        size_t current_autocor_offset = 1;
        for (size_t i = 0; i != 3; ++i)
        {
            //sample more and more thinly, starting from every 1'th
            this->sampler->sample(this->tuning_num_points, 0, cumul_autocor_offset,
                                  (verbose ? & proposal_mean : NULL), // only calculate extra stats if in verbose mode
                                  (verbose ? & proposal_variance : NULL),
                                  sample_points_buf);
            
            current_autocor_offset =
                best_autocorrelation_offset(sample_points_buf,
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
                        locus->position);
                fflush(stdout);
            }
        }

        //here it doesn't make sense to take all of the samples if
        //the best_autocor offset isn't 1
        multivariate_mean(sample_points_buf, full_ndim,
                          this->tuning_num_points, estimated_mean);

        if (verbose)
        {
            fprintf(stdout, "MH: %Zu, position: %i, proposal mean: %g, "
                    "proposal variance: %g", cumul_autocor_offset,
                    locus->position, proposal_mean, proposal_variance);
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
                                                   this->model->composition_prior_alphas);
        
    }
    return cumul_autocor_offset;
}

// uses sample_points_buf as a scratch space for tuning the slice_sampler
size_t posterior_wrapper::tune_ss(double *sample_points_buf)
{
    //SliceSampling slice_sampler(truncated_ndim, num_bits_per_dim, is_log_integrand, 1);

    this->slice_sampler->Initialize();
    this->slice_sampler->sample(this->model, 
                                this->mode_point,
                                initial_sampling_range, 1, sample_points_buf,
                                tuning_num_points);
    
    const size_t truncated_ndim = 3;

    size_t best_autocor_offset =
        best_autocorrelation_offset(sample_points_buf, truncated_ndim, 
                                    tuning_num_points, autocor_max_offset, 
                                    autocor_max_value);
    return best_autocor_offset;
}


void posterior_wrapper::tune(sample_details *sd)
{
    // try MH first.
    sd->autocor_offset = this->tune_mh(sd->locus, sd->sample_points);

    if (sd->autocor_offset <= this->target_autocor_offset)
        sd->samp_method = METROPOLIS_HASTINGS;
    else
    {
        sd->autocor_offset = this->tune_ss(sd->sample_points);
        sd->samp_method = (sd->autocor_offset <= this->target_autocor_offset)
            ? SLICE_SAMPLING
            : FAILED;
    }
}


void posterior_wrapper::sample(sample_details *sd, size_t num_points)
{
    switch(sd->samp_method)
    {
    case METROPOLIS_HASTINGS:
        this->sampler->sample(num_points, 0, sd->autocor_offset, NULL, NULL, sd->sample_points);
        break;
    case SLICE_SAMPLING:
        this->slice_sampler->sample(this->model, 
                                    this->mode_point,
                                    this->initial_sampling_range, 
                                    sd->autocor_offset,
                                    sd->sample_points,
                                    num_points);
        break;
    case FAILED:
        break;
    }
}

// generate sample points according to the parameters set in posterior_wrapper
// alt_sample_points must be of size this->final_num_points * 4
// write 
void posterior_wrapper::sample(PileupSummary *locus, double *sample_points_buf, char *algorithm_used)
{
    size_t cumul_autocor_offset = this->tune_mh(locus, sample_points_buf);
    // double proposal_mean, proposal_variance;

    if (cumul_autocor_offset <= this->target_autocor_offset)
    {
        //metropolis hastings succeeded
        this->sampler->sample(this->final_num_points, 0, cumul_autocor_offset, NULL, NULL, sample_points_buf);
        strcpy(algorithm_used, "MH");
    }
    else
    {
        //metropolis hastings failed.  try slice sampling
        size_t best_autocor_offset = this->tune_ss(sample_points_buf);
        
        if (best_autocor_offset <= this->target_autocor_offset)
        {
            //slice sampling succeeded
            this->slice_sampler->sample(this->model, 
                                        this->mode_point,
                                        this->initial_sampling_range, 
                                        best_autocor_offset,
                                        sample_points_buf,
                                        this->final_num_points);
            strcpy(algorithm_used, "SS");
        }
        else
        {
            strcpy(algorithm_used, "--");
        }
    }
}


// generate a set of normalized values of the posterior at the
// specified points
void posterior_wrapper::values(double *points, size_t num_points,
                               double *values)
{
    double *point = points;
    double *end = points + (num_points * 4);
    double *val = values;
    double *vend = values + num_points;

    for (; point != end; point += 4, ++val)
    {
        *val = this->model->log_likelihood(point);
        assert(! isinf(*val));
        assert(! isnan(*val));
    }
    // transform these into normalized values
    double maxv = values[0];
    for (val = values; val != vend; ++val)
    {
        maxv = std::max(maxv, *val);
    }
    double sum = 0.0;
    for (val = values; val != vend; ++val)
    {
        gsl_sf_result res;
        int err = gsl_sf_exp_e(*val - maxv, &res);
        if (err != 0 && err != GSL_EUNDRFLW)
        {
            fprintf(stderr, "Error: This should either be no error or an underflow\n");
            exit(10);
        }

        *val = res.val;
        assert(! isnan(*val));
        sum += *val;
    }
    for (val = values; val != vend; ++val)
    {
        *val /= sum;
    }    
    
    
}


char *posterior_wrapper::print_quantiles(sample_details *sd, char *out_buffer)
{
    char line_label[2048];
    size_t effective_depth = sd->locus->read_depth; // !!! fix this

    sprintf(line_label,
            "%s\t%c\t%s\t%i\t%c\t%Zu\t%Zu",
            this->label_string, 
            (char)sd->samp_method, 
            sd->locus->reference, 
            sd->locus->position, 
            sd->locus->reference_base, 
            sd->locus->read_depth, 
            effective_depth
            );

    const char *dimension_labels[] = { "A", "C", "G", "T" };

    out_buffer = 
        print_marginal_quantiles(out_buffer,
                                 sd->sample_points,
                                 this->final_num_points,
                                 this->mode_point,
                                 line_label, 
                                 dimension_labels, 
                                 "+", 
                                 this->quantiles,
                                 this->num_quantiles);
    return out_buffer;
}


// this should correspond to print_discrete_comp
size_t discrete_comp_locus_bytes(size_t num_discrete_values)
{
    return 30 + 10 + 1 + 10 + 4 + 1 + ((4 + 1 + 7) * num_discrete_values);
}


struct first_comp_more
{
    bool operator()(std::pair<double, size_t> const& a,
                    std::pair<double, size_t> const& b)
    {
        return b.first < a.first;
    }
};

// prints an abbreviated format of base composition as represented by
// a discrete set of posterior points.
char *print_discrete_comp(PileupSummary *locus,
                          const char *sample_label,
                          double *discrete_values,
                          size_t num_discrete_values,
                          size_t num_items_to_print,
                          double min_value_to_print,
                          char *out_buf)
{

    std::pair<double, size_t> *val_ord = new std::pair<double, size_t>[num_discrete_values];
    std::pair<double, size_t> *vo = val_ord;
    double *val = discrete_values;
    for (size_t d = 0; d != num_discrete_values; ++d)
    {
        (*vo).first = *val;
        (*vo).second = d;
        ++vo;
        ++val;
    }
    first_comp_more fcm;
    std::partial_sort(val_ord, val_ord + num_items_to_print, val_ord + num_discrete_values, fcm);

    out_buf += sprintf(out_buf, 
                       "%s\t%s\t%i\t%c\t%Zu\t", 
                       sample_label,
                       locus->reference,
                       locus->position,
                       locus->reference_base,
                       locus->read_depth);
    
    const char *sep = "";
    vo = val_ord;
    for (size_t d = 0; d != num_items_to_print; ++d)
    {
        if ((*vo).first > min_value_to_print)
        {
            out_buf += sprintf(out_buf, "%s%Zu:%-4.3f", sep, (*vo).second, (*vo).first);
            sep = ",";
        }
        ++vo;
    }
    out_buf += sprintf(out_buf, "\n");

    delete val_ord;

    return out_buf;
}


// returns next write position after writing to out_buffer
// sample_points_buf is a scratch space for writing sample points
// must be this->final_num_points * 4 size
char *posterior_wrapper::process_line_comp(const char *pileup_line,
                                           char *out_buffer,
                                           double *sample_points_buf,
                                           float test_quantile,
                                           float min_test_quantile_value)
{

    PileupSummary locus;
    locus.load_line(pileup_line);
    locus.parse(this->min_quality_score);
    this->model->locus_data = &locus.counts;
    this->params->pack(&locus.counts);
    this->find_mode();

    // first-pass test.  if the mode point indicates that there is no
    // chance of a higher confidence non-reference base composition, 
    // skip all further tests and do not output anything.
    const char *nucs = "ACGTacgt", *query;
    int ref_ind = (query = strchr(nucs, locus.reference_base)) ? (query - nucs) % 4 : -1;
    if (ref_ind >= 0 && 
        this->mode_point[ref_ind] < 1.0 - min_test_quantile_value)
        return out_buffer; 

    sample_details sd;
    sd.locus = &locus;
    sd.sample_points = sample_points_buf;
    this->tune(&sd);
    this->sample(&sd, this->final_num_points);

    // second-pass test.
    double d_test_quantile = (double)test_quantile;
    double test_quantile_value;
    if (ref_ind >= 0)
        compute_marginal_quantiles(sd.sample_points,
                                   this->final_num_points,
                                   ref_ind,
                                   &d_test_quantile,
                                   1,
                                   &test_quantile_value);

    
    // invert the value threshold.  we want to test for distance of
    // reference base away from 1, rather than presence of
    // non-reference base.
    if (ref_ind >= 0 && test_quantile_value >= (1.0 - min_test_quantile_value))
        return out_buffer;


    // if the reference base is not in "ACGTacgt", we always print the
    // locus.
    out_buffer = this->print_quantiles(&sd, out_buffer);
    
    if (cdfs_output_fh != NULL)
    {
        const size_t full_ndim = 4;
        pthread_mutex_lock(file_writing_mutex);
        print_numerical_cdfs(cdfs_output_fh, 
                             this->label_string, 
                             sample_points_buf,
                             this->final_num_points, 
                             full_ndim);
        
        pthread_mutex_unlock(file_writing_mutex);
    }
    
    return out_buffer;
}



char *posterior_wrapper::process_line_mode(const char *pileup_line,
                                           char *out_buffer)
{

    PileupSummary locus;
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
    this->model->locus_data = &locus.counts;
    // this->posterior->initialize(this->gradient_tolerance, this->max_modefinding_iterations, 
    //                             initial_point, verbose);

    this->model->find_mode_point(this->gradient_tolerance,
                                 this->max_modefinding_iterations,
                                 this->initial_point,
                                 this->on_zero_boundary,
                                 this->verbose,
                                 this->mode_point);

    out_buffer += 
        sprintf(out_buffer, 
                "%s\t%s\t%i\t%c\t%Zu\t%Zu\t%5.5f\t%5.5f\t%5.5f\t%5.5f\n",
                label_string, 
                locus.reference, 
                locus.position, 
                locus.reference_base, 
                locus.read_depth,
                effective_depth,
                this->mode_point[0],
                this->mode_point[1],
                this->mode_point[2],
                this->mode_point[3]
                );
    
    return out_buffer;
}




void *comp_worker(void *args)
{
    wrapper_input *input = static_cast<wrapper_input *>(args);
    std::vector<char *>::iterator it;
    char **out = input->out_start;
    double *sample_points_buf = new double[input->worker->final_num_points * 4];
    for (it = input->beg; it != input->end; ++it)
    {
        input->worker->process_line_comp(*it, *out, sample_points_buf,
                                         input->test_quantile,
                                         input->min_test_quantile_value);
        ++out;
    }
    delete sample_points_buf;
    pthread_exit((void*) 0);
}



void *mode_worker(void *args)
{
    wrapper_input *input = static_cast<wrapper_input *>(args);
    std::vector<char *>::iterator it;
    char **out = input->out_start;
    for (it = input->beg; it != input->end; ++it)
    {
        input->worker->process_line_mode(*it, *out);
        ++out;
    }
    pthread_exit((void*) 0);
}
