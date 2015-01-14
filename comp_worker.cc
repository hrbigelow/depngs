#include "comp_worker.h"

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
#include "defs.h"
#include "cache.h"
#include "locus.h"

#if 0
posterior_wrapper::posterior_wrapper(const char *jpd_data_params_file,
                                     double *prior_alphas,
                                     size_t min_quality_score,
                                     double *quantiles,
                                     size_t n_quantiles,
                                     const char *label_string,
                                     FILE *cdfs_output_fh,
                                     pthread_mutex_t *file_writing_mutex,
                                     bool verbose) :
    s(s),
    verbose(verbose),
    min_quality_score(min_quality_score),
    n_quantiles(n_quantiles),
    cdfs_output_fh(cdfs_output_fh),
    file_writing_mutex(file_writing_mutex)
{
    this->label_string = strdup(label_string);

    this->quantiles = new double[this->n_quantiles];
    std::copy(quantiles, quantiles + n_quantiles, this->quantiles);

    int i;
    for (i = 0; i != 4; ++i)
        this->initial_point[i] = 0.25;

    //this->prior_alpha0 = std::accumulate(prior_alphas, prior_alphas + 4, 0.0);
    
    this->params = new NucleotideStats();
    this->params->initialize(jpd_data_params_file);
    this->model = new ErrorEstimate();
    this->model->set_prior_alphas(prior_alphas);
    this->model->model_params = this->params;
    // this->prior = new Dirichlet();
    this->sampler = new Metropolis(this->model, new Dirichlet(),
                                   this->s.final_n_points);
    this->slice_sampler = new SliceSampling(1);
}



posterior_wrapper::~posterior_wrapper()
{
    free(this->label_string);
    delete this->quantiles;
    delete this->sampler;
    delete this->slice_sampler;
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
        this->model->find_mode_point(this->s.gradient_tolerance,
                                     this->s.max_modefinding_iterations,
                                     this->initial_point,
                                     this->on_zero_boundary,
                                     this->verbose,
                                     this->mode_point);
}


// estimate the alphas of a dirichlet that mimics the posterior based
// on raw counts.
void alphas_from_raw_counts(packed_counts *counts, double *prior_alpha, double *est_alpha)
{
    char base;
    size_t qual, strand;
    memcpy(est_alpha, prior_alpha, sizeof(double) * NUM_NUCS);
    for (size_t d = 0; d != counts->num_data; ++d)
    {
        Nucleotide::decode(counts->stats_index[d], &base, &qual, &strand);
        est_alpha[Nucleotide::base_to_index[(int)base]] += counts->raw_counts[d];
    }
}


// tune the posterior, return the cumulative autocorrelation offset
// sample_points_buf is used as a space for writing sample points
// during the tuning process.
size_t posterior_wrapper::tune_mh(PileupSummary *locus, double *sample_points_buf,
                                  double *estimated_mean /* output */)
{

    // start out by setting the proposal alphas to correspond with the
    // counts of basecalls in the locus plus the prior alphas.
    // !!! this is inefficient...
    double alpha[NUM_NUCS];

    alphas_from_raw_counts(&locus->counts, this->model->prior_alpha, alpha);

    this->sampler->proposal->update(alpha);
    
    // set the current point to the normalized alphas
    for (size_t n = 0; n != NUM_NUCS; ++n)
        estimated_mean[n] = this->sampler->current_point[n] = 
            this->sampler->proposal->alpha[n] /
            this->sampler->proposal->alpha0;

    
    // set the proposal alphas
    this->sampler->proposal->set_alphas_from_mean_or_bound(this->sampler->current_point,
                                                           this->model->prior_alpha);

    //metropolis hastings
    double proposal_mean, proposal_variance;
    size_t cumul_autocor_offset = 0;

    
    for (size_t iter = 0; iter != this->s.max_tuning_iterations; ++iter)
    {
        cumul_autocor_offset = 1;
        size_t current_autocor_offset = 1;
        for (size_t i = 0; i != 3; ++i)
        {
            // start the sampling procedure from the mean (doesn't really matter too much)
            //sample more and more thinly, starting from every 1'th
            this->sampler->sample(this->s.tuning_n_points, 0, cumul_autocor_offset,
                                  estimated_mean,
                                  (verbose ? & proposal_mean : NULL), // only calculate extra stats if in verbose mode
                                  (verbose ? & proposal_variance : NULL),
                                  sample_points_buf);
            
            current_autocor_offset =
                best_autocorrelation_offset(sample_points_buf,
                                            NUM_NUCS,
                                            this->s.tuning_n_points,
                                            this->s.autocor_max_offset, 
                                            this->s.autocor_max_value);

            if (current_autocor_offset == 1)
                break;

            cumul_autocor_offset *= current_autocor_offset;
            if (verbose)
            {
                fprintf(stdout, "MH: current: %Zu, cumulative: %Zu, position: %i\n", 
                        current_autocor_offset, cumul_autocor_offset, 
                        locus->position);
                fflush(stdout);
            }
        }

        if (verbose)
        {
            fprintf(stdout, "MH: %Zu, position: %i, proposal mean: %g, "
                    "proposal variance: %g", cumul_autocor_offset,
                    locus->position, proposal_mean, proposal_variance);
            fprintf(stdout, ", estimated mean:");

            multivariate_mean(sample_points_buf, NUM_NUCS,
                              this->s.tuning_n_points, estimated_mean);

            for (size_t d = 0; d != NUM_NUCS; ++d)
            {
                fprintf(stdout, "\t%g", estimated_mean[d]);
            }
            fprintf(stdout, "\n");
            fflush(stdout);
        }

        if (cumul_autocor_offset <= this->s.target_autocor_offset)
            break;

        // ideally, i had wanted to just use this, but for some
        // reason, the estimated mean can drift into having alphas
        // close to zero.  theoretically, this shouldn't happen
        // because the posterior comprises additional observations
        // besides the prior.

        // this->sampler->proposal->set_alphas_from_mean(estimated_mean);
        //here it doesn't make sense to take all of the samples if
        //the best_autocor offset isn't 1
        multivariate_mean(sample_points_buf, NUM_NUCS,
                          this->s.tuning_n_points, estimated_mean);

        this->sampler->proposal->set_alphas_from_mean_or_bound(estimated_mean,
                                                               this->model->prior_alpha);
        
    }
    return cumul_autocor_offset;
}

// uses sample_points_buf as a scratch space for tuning the slice_sampler
size_t posterior_wrapper::tune_ss(PileupSummary *locus, 
                                  double *sample_points_buf,
                                  double *estimated_mean)
{
    //SliceSampling slice_sampler(truncated_ndim, num_bits_per_dim, is_log_integrand, 1);

    this->slice_sampler->Initialize();

    double alpha0 = 0;
    alphas_from_raw_counts(&locus->counts, this->model->prior_alpha, estimated_mean);

    for (size_t n = 0; n != NUM_NUCS; ++n)
        alpha0 += estimated_mean[n];

    for (size_t n = 0; n != NUM_NUCS; ++n)
        estimated_mean[n] /= alpha0;

    this->slice_sampler->sample(this->model, 
                                estimated_mean,
                                this->s.initial_sampling_range, 1, sample_points_buf,
                                this->s.tuning_n_points);
    
    const size_t truncated_ndim = 3;

    size_t best_autocor_offset =
        best_autocorrelation_offset(sample_points_buf, truncated_ndim, 
                                    this->s.tuning_n_points, 
                                    this->s.autocor_max_offset, 
                                    this->s.autocor_max_value);
    return best_autocor_offset;
}


void posterior_wrapper::tune(sample_details *sd, double *estimated_mean)
{
    // try MH first.
    sd->autocor_offset = this->tune_mh(&sd->locus, sd->sample_points, estimated_mean);

    if (sd->autocor_offset <= this->s.target_autocor_offset)
        sd->samp_method = METROPOLIS_HASTINGS;
    else
    {
        sd->autocor_offset = this->tune_ss(&sd->locus, sd->sample_points, estimated_mean);
        sd->samp_method = (sd->autocor_offset <= this->s.target_autocor_offset)
            ? SLICE_SAMPLING
            : FAILED;
    }
}


void posterior_wrapper::sample(sample_details *sd, double *initial_point, size_t n_points)
{
    // burn_in is usually only necessary with 
    size_t burn_in = 10;
    switch(sd->samp_method)
    {
    case METROPOLIS_HASTINGS:
        this->sampler->sample(n_points, burn_in, sd->autocor_offset,
                              initial_point,
                              NULL, NULL, sd->sample_points);
        break;
    case SLICE_SAMPLING:
        this->slice_sampler->sample(this->model, 
                                    initial_point,
                                    this->s.initial_sampling_range, 
                                    sd->autocor_offset,
                                    sd->sample_points,
                                    n_points);
        break;
    case FAILED:
        break;
    }
}
#endif


// generate a set of normalized values of the posterior at the
// specified points
#if 0
void posterior_wrapper::values(double *points, size_t n_points,
                               double *values)
{
    double *point = points;
    double *end = points + (n_points * 4);
    double *val = values;
    double *vend = values + n_points;

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
#endif


char *print_quantiles(const double *quantiles,
                      size_t n_quantiles,
                      const char *label_string,
                      struct locus_sampling *ls, 
                      char *out_buffer)
{
    char line_label[2048];

    sprintf(line_label,
            "%s\t%s\t%i\t%c\t%Zu\t%Zu",
            label_string, 
            ls->locus.reference, 
            ls->locus.position, 
            ls->locus.reference_base, 
            ls->locus.read_depth, 
            ls->locus.read_depth_high_qual
            );

    out_buffer = print_marginal_quantiles(out_buffer,
                                          ls->sample_points,
                                          ls->n_sample_points,
                                          line_label,
                                          "+", 
                                          quantiles,
                                          n_quantiles);
    return out_buffer;
}



// prints an abbreviated format of base composition as represented by
// a discrete set of posterior points.
#if 0
char *print_discrete_comp(PileupSummary *locus,
                          const char *sample_label,
                          double *discrete_values,
                          size_t n_discrete_values,
                          size_t n_items_to_print,
                          double min_value_to_print,
                          char *out_buf)
{

    std::pair<double, size_t> *val_ord = new std::pair<double, size_t>[n_discrete_values];
    std::pair<double, size_t> *vo = val_ord;
    double *val = discrete_values;
    for (size_t d = 0; d != n_discrete_values; ++d)
    {
        (*vo).first = *val;
        (*vo).second = d;
        ++vo;
        ++val;
    }
    first_comp_more fcm;
    std::partial_sort(val_ord, val_ord + n_items_to_print, val_ord + n_discrete_values, fcm);

    out_buf += sprintf(out_buf, 
                       "%s\t%s\t%i\t%c\t%Zu\t", 
                       sample_label,
                       locus->reference,
                       locus->position,
                       locus->reference_base,
                       locus->read_depth);
    
    const char *sep = "";
    vo = val_ord;
    for (size_t d = 0; d != n_items_to_print; ++d)
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
#endif

struct sample_attributes
make_sample_attributes(const char *jpd_file,
                       const char *label_string,
                       const char *pileup_file)
{
    struct sample_attributes s;
    s.nuc_stats = make_nucleotide_stats();
    nucleotide_stats_initialize(jpd_file, &s.nuc_stats);
    if (strlen(label_string) + 1 > sizeof(s.label_string)) 
    {
        fprintf(stderr, "%s: error: sample label string must be "
                "less than %Zu characters\n", __func__,
                sizeof(s.label_string));
        exit(1);
    }
    strcpy(s.label_string, label_string);
    s.fh = fopen(pileup_file, "r");
    if (! s.fh)
    {
        fprintf(stderr, "%s: error: couldn't open pileup input file %s\n",
                __func__, pileup_file);
        exit(1);
    }
    return s;
}


/* initialize the locus in 'ls' defined by sd->current */
inline void init_pileup_locus(const struct nucleotide_stats *stats,
                              size_t min_quality_score,
                              locus_sampling *ls)
{
    ls->locus.load_line(ls->current);
    ls->locus_ord = init_locus(ls->current);
    ls->locus.parse(min_quality_score);
    ls->n_sample_points = 0;
    // input->worker[s]->model->locus_data = &ls->locus.counts;
    nucleotide_stats_pack(stats, &ls->locus.counts);
}


/* advance ls->current, then initialize if there is another locus */
void refresh_locus(const struct nucleotide_stats *stats,
                   size_t min_quality_score,
                   locus_sampling *ls)
{
    assert(ls->current != ls->end);
    if ((ls->current = strchr(ls->current, '\n') + 1) == ls->end) 
        ls->is_next = false;

    else
    {
        init_pileup_locus(stats, min_quality_score, ls);
        ls->dist_printed = 0;
    }
}


/* returns next write position after writing to out_buffer
   sample_points_buf is a scratch space for writing sample points must
   be this->final_n_points * 4 size */
char *process_line_comp(struct comp_worker_input *cw,
                        struct locus_sampling *ls,
                        char *out_buffer,
                        float test_quantile,
                        float min_test_quantile_value)
{
    const struct cpd_count
        *pc = ls->locus.counts.stats,
        *pce = pc + ls->locus.counts.num_data;
    while (pc != pce)
    {
        fprintf(stderr, "{ %g, %g, %g, %g, %lu },\n",
        pc->cpd[0], pc->cpd[1], pc->cpd[2], pc->cpd[3], pc->ct);
        ++pc;
    }

    double 
        test_quantile_value,
        ref_test_quantile = 1.0 - (double)test_quantile,
        max_test_quantile_value = 1.0 - min_test_quantile_value;

    // first-pass test.  if reference-base component of mode point is
    // too close to 1, then the test quantile value will be as well.
    char nucs[] = "ACGTacgt", *query;
    int ref_ind = (query = strchr(nucs, ls->locus.reference_base)) ? (query - nucs) % 4 : -1;

#ifndef _SKIP_FIRST_PASS_TEST_
    // if (ref_ind >= 0 && 
    //     this->mode_point[ref_ind] > max_test_quantile_value)
    //     return out_buffer; 
#endif

    /* opportunity here to not re-sample points if the tuning has okay
       autocorrelation */
    double estimated_mean[NUM_NUCS], proposal_alpha[NUM_NUCS];;
    size_t cumul_aoff = 
        tune_proposal(&ls->locus.counts,
                      &cw->pset, proposal_alpha, estimated_mean,
                      ls->sample_points);

    // this->tune(&ls, initial_point);
    metropolis_sampling(0, cw->pset.final_n_points, &ls->locus.counts,
                        cw->pset.logu, proposal_alpha, cumul_aoff,
                        ls->sample_points);
    // this->sample(&ls, initial_point, this->s.final_n_points);

    // we always print a base composition estimate for loci with
    // ref=N.  So, we only do the filtering test if ref!=N
    if (ref_ind >= 0)
        compute_marginal_quantiles(ls->sample_points,
                                   cw->pset.final_n_points,
                                   ref_ind,
                                   &ref_test_quantile,
                                   1,
                                   &test_quantile_value);

    
    // invert the value threshold.  we want to test for distance of
    // reference base away from 1, rather than presence of
    // non-reference base.
    if (ref_ind >= 0 && test_quantile_value > max_test_quantile_value)
        return out_buffer;


    // if the reference base is not in "ACGTacgt", we always print the
    // locus.
    out_buffer = print_quantiles(cw->quantiles,
                                 cw->n_quantiles,
                                 cw->sample_atts.label_string, 
                                 ls, out_buffer);
    return out_buffer;
}


/* iterate through the input range of lines.  nullifies the newline
   characters of the input.  responsible for allocating output
   buffer. conforms to thread_queue_worker_t.  in_buf and out_buf
   point to arrays of size 1 in this case. */
void comp_worker(void *par, 
                 const struct managed_buf *in_buf,
                 struct managed_buf *out_buf)
{
    struct comp_worker_input *cw = 
        (struct comp_worker_input *)par;

    char *write_ptr = out_buf->buf;
    size_t max_line = 1000;
    struct locus_sampling ls = {
        PileupSummary(),
        false,
        0,
        (double *)malloc(cw->pset.final_n_points * sizeof(double)),
        0,
        0,
        in_buf->buf,
        in_buf->buf + in_buf->size,
        min_pair_ord
    };
    
    init_pileup_locus(&cw->sample_atts.nuc_stats,
                      cw->pset.min_quality_score,
                      &ls);

    while (ls.current != ls.end)
    {
        ALLOC_GROW_REMAP(out_buf->buf, write_ptr, 
                         write_ptr - out_buf->buf + max_line, out_buf->alloc);
        write_ptr =
            process_line_comp(cw, &ls, write_ptr,
                              cw->test_quantile,
                              cw->min_test_quantile_value);

        refresh_locus(&cw->sample_atts.nuc_stats,
                      cw->pset.min_quality_score, &ls);
    }
    out_buf->size = write_ptr - out_buf->buf;
    free(ls.sample_points);
}
