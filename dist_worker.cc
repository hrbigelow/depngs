#include "dist_worker.h"
#include "pileup_tools.h"
#include "comp_functor.h"
#include "locus_comp.h"
#include "error_estimate.h"
#include "vcf.h"

#include <cstdio>
#include <cstring>
#include <cassert>
#include <gsl/gsl_math.h>
#include <algorithm>


eval_dist_matrix::eval_dist_matrix() :
    points(NULL),
    num_points(0),
    dist_index(NULL),
    num_distances(0),
    distances(NULL),
    inclusion_threshold(0),
    max_points_to_print(0),
    min_value_to_print(0)
{}


dist_worker_input::dist_worker_input(size_t thread_num,
                                     size_t num_samples,
                                     size_t num_sample_pairs,
                                     size_t num_sample_point_pairings,
                                     double *dist_quantiles,
                                     size_t num_dist_quantiles,
                                     double *comp_quantiles,
                                     size_t num_comp_quantiles,
                                     eval_dist_matrix *lattice,
                                     char *out_dist,
                                     char *out_comp,
                                     char *out_vcf,
                                     size_t *pair_sample1,
                                     size_t *pair_sample2,
                                     std::map<const char *, size_t, ltstr> *contig_order) :
    thread_num(thread_num),
    num_samples(num_samples), num_sample_pairs(num_sample_pairs),
    num_sample_point_pairings(num_sample_point_pairings),
    dist_quantiles(dist_quantiles), num_dist_quantiles(num_dist_quantiles),
    comp_quantiles(comp_quantiles), num_comp_quantiles(num_comp_quantiles),
    lattice(lattice),
    out_dist(out_dist), out_comp(out_comp),
    out_vcf(out_vcf),
    pair_sample1(pair_sample1),
    pair_sample2(pair_sample2), contig_order(contig_order)
{ 
    this->worker = new posterior_wrapper *[this->num_samples];
    this->beg = new std::vector<char *>::iterator[this->num_samples];
    this->end = new std::vector<char *>::iterator[this->num_samples];
    this->less_locus.contig_order = this->contig_order;
    this->equal_locus.contig_order = this->contig_order;
    this->randgen = gsl_rng_alloc(gsl_rng_taus);
    gsl_rng_set(this->randgen, 0);
}


dist_worker_input::dist_worker_input() :
    worker(NULL), num_samples(0), num_sample_pairs(0),
    num_sample_point_pairings(0),
    dist_quantiles(NULL), num_dist_quantiles(0),
    comp_quantiles(NULL), num_comp_quantiles(0),
    beg(NULL),
    end(NULL),
    out_dist(NULL), out_comp(NULL),
    out_discomp(NULL),
    out_vcf(NULL),
    pair_sample1(NULL),
    pair_sample2(NULL),
    contig_order(NULL)
{
}


dist_worker_input::~dist_worker_input()
{
    for (size_t s = 0; s != this->num_samples; ++s)
    {
        if (this->worker[s] != NULL)
        {
            delete this->worker[s];
        }
    }
    if (this->worker != NULL) delete this->worker;
    if (this->beg != NULL) delete this->beg;
    if (this->end != NULL) delete this->end;
    gsl_rng_free(this->randgen);
}

// assumes the chromosome has a number


size_t distance_quantiles_locus_bytes(size_t num_quantiles)
{
    return 3 + 3 + 10 + 10 + (7 * num_quantiles) + (4 + num_quantiles);
}



// compute and print out distance quantiles, based on quasi-random
// pairings of two samplings for efficiency, the 'random' pairing is
// done simply by cycling through both sets of sample points, but
// starting in the middle for the second set.
// return the next write position.


char *print_distance_quantiles(const double *points1,
                               const double *points2,
                               const char *contig,
                               size_t position,
                               dist_worker_input *wi,
                               char *out_dist_buf,
                               size_t pair_index)
{
    size_t s1 = wi->pair_sample1[pair_index];
    size_t s2 = wi->pair_sample2[pair_index];

    // if (position == 58526522)
    // {
    //     size_t x = 0;
    //     ++x;
    // }

    size_t num_random_pairs = wi->num_sample_point_pairings;
    size_t num_sample_points1 = wi->worker[s1]->final_num_points;
    size_t num_sample_points2 = wi->worker[s2]->final_num_points;
    double *square_dist = new double[num_random_pairs];

    size_t pi1, pi2 = 0;
    for (size_t p = 0; p != num_random_pairs; ++p)
    {
        pi1 = (p % num_sample_points1) * 4;
        pi2 = gsl_rng_uniform_int(wi->randgen, num_sample_points2) * 4;
        square_dist[p] = 
            gsl_pow_2(points1[pi1] - points2[pi2])
            + gsl_pow_2(points1[pi1 + 1] - points2[pi2 + 1])
            + gsl_pow_2(points1[pi1 + 2] - points2[pi2 + 2])
            + gsl_pow_2(points1[pi1 + 3] - points2[pi2 + 3]);

        // fprintf(stderr, 
        //         "%s\t%s\t%s\t%Zu\t%4.3f\n",
        //         wi->worker[s1]->label_string,
        //         wi->worker[s2]->label_string,
        //         contig, position,
        //         sqrt(square_dist[p]));
                
    }

    // find the quantiles
    double *start = square_dist;
    double *end = square_dist + num_random_pairs;
    double *cut;
    double qval;

    out_dist_buf += sprintf(out_dist_buf, "%s\t%s\t%s\t%Zu", 
                            wi->worker[s1]->label_string,
                            wi->worker[s2]->label_string,
                            contig, position);
    
    for (size_t q = 0; q != wi->num_dist_quantiles; ++q)
    {
        cut = square_dist + static_cast<size_t>(std::round(wi->dist_quantiles[q] * num_random_pairs));
        std::nth_element(start, cut, end);
        qval = (cut == end) ? -1.0 : sqrt(*cut);
        start = cut;
        out_dist_buf += sprintf(out_dist_buf, "\t%7.4f", qval);
    }

    out_dist_buf += sprintf(out_dist_buf, "\n");

    delete square_dist;

    return out_dist_buf;
}




// print out distance quantiles for the next locus, for all pairs
char *next_distance_quantiles_aux(dist_worker_input *input, 
                                  sample_details *sd,
                                  size_t global_s,
                                  const double *null_points,
                                  char *out_buf)
{
    if (out_buf == NULL)
    {
        return out_buf;
    }

    char contig[100];
    size_t position;
    double *mode1, *mode2;

    sscanf(*sd[global_s].current, "%s\t%zu", contig, &position);

    for (size_t p = 0; p != input->num_sample_pairs; ++p)
    {
        size_t s1 = input->pair_sample1[p];
        size_t s2 = input->pair_sample2[p];
        sample_details &sd1 = sd[input->pair_sample1[p]];
        sample_details &sd2 = sd[input->pair_sample2[p]];

        mode1 = input->worker[s1]->mode_point, end1 = mode1 + 4;
        mode2 = input->worker[s2]->mode_point;

        float mode_square_dist = 0;
        while (mode1 != end)
            mode_square_dist += *mode1++ * *mode2++;

        if (mode_square_dist < input->min_high_conf_dist_squared)
        {
            // modes don't differ enough
            continue;
        }
        // preliminary sampling
        if (sd1.is_next && ! sd1.num_sample_points)
        {
            input->worker[s1]->tune(sd1.locus, sd1.sample_points,
                                    &sd1.samp_method, &sd1.autocor_offset);

            input->worker[s1]->sample(sd1.locus, 
                                      sd1.sample_points, 
                                      input->num_prelim_points,
                                      sd1.autocor_offset,
                                      sd1.samp_method);
            sd1.num_sample_points = input->num_prelim_points;
        }
        if (sd2.is_next && ! sd2.num_sample_points)
        {
            input->worker[s2]->tune(sd2.locus, sd2.sample_points,
                                    &sd2.samp_method, &sd2.autocor_offset);

            input->worker[s2]->sample(sd2.locus, 
                                      sd2.sample_points, 
                                      input->num_prelim_points,
                                      sd2.autocor_offset,
                                      sd2.samp_method);
            sd2.num_sample_points = input->num_prelim_points;
        }

        

        const double *points1 = sd1.is_next ? sd1.sample_points : null_points;
        const double *points2 = sd2.is_next ? sd2.sample_points : null_points;
        
        out_buf = 
            print_distance_quantiles(points1, points2, contig, position, input, out_buf, p);
        }
    }
    return out_buf;
}


// find global_s, and init the is_next field for all samples
// this is run
size_t init_global_and_next(dist_worker_input *input, sample_details *samples)
{
    size_t global_s = 0;
    for (size_t s = 0; s != input->num_samples; ++s)
    {
        if (samples[s].current != input->end[s]
            && (samples[global_s].current == input->end[global_s]
                || input->less_locus(*samples[s].current, *samples[global_s].current)
                ))
        {
            global_s = s;
        }
    }
    bool global_at_end = (samples[global_s].current == input->end[global_s]);
    for (size_t s = 0; s != input->num_samples; ++s)
    {
        samples[s].is_next = 
            (samples[s].current != input->end[s]) 
            && (! global_at_end)
            && input->equal_locus(*samples[s].current, *samples[global_s].current);
    }

    return global_s;
}


// initialize the locus defined by sd->current.
// assume sd->current is a valid iterator
void init_locus(dist_worker_input *input,
                size_t sample_id,
                sample_details *sd)
{
    size_t s = sample_id;
    sd->locus = new PileupSummary(0);
    sd->locus->load_line(*sd->current);
    sd->locus->parse(input->worker[s]->min_quality_score);

    input->worker[s]->model->locus_data = &sd->locus->counts;
    input->worker[s]->params->pack(&sd->locus->counts);
                    
    this->worker[s]->find_mode();
}


// optionally destroy previous locus associated with sd,
// and create a new one, initializing it with sample points
// assume sd->locus and sd->current are initialized
void refresh_locus(dist_worker_input *input,
                   size_t sample_id,
                   sample_details *sd)
{
    size_t s = sample_id;
    assert(sd->current != input->end[s]);

    if (sd->locus != NULL)
    {
        delete sd->locus;
        sd->locus = NULL;
    }
    if (++(sd->current) == input->end[s])
    {
        sd->locus = NULL;
    }
    else
    {
        init_locus(input, sample_id, sd);
    }
}


// refresh any loci that are marked as 'next'
// does not initialize is_next for new loci.  (that can only be done
void advance_loci_aux(dist_worker_input *input,
                      sample_details *sd,
                      char **out_comp_buf,
                      char **out_vcf_buf)
{

    if (*out_vcf_buf != NULL)
    {
        *out_vcf_buf = print_vcf_line(sd, input, *out_vcf_buf);
    }

    for (size_t s = 0; s != input->num_samples; ++s)
    {
        if (! sd[s].is_next)
        {
            continue;
        }

        if (*out_comp_buf != NULL && sd[s].num_sample_points)
        {
            *out_comp_buf = 
                input->worker[s]->print_quantiles(sd[s].locus, sd[s].algorithm_used,
                                                  sd[s].sample_points, *out_comp_buf);


            /*
            for (size_t p = 0; p != input->worker[s]->final_num_points; ++p)
            {
                // sample_id, chromosome, position, fA, fC, fG, fT
                fprintf(stdout, "%s\t%s\t%i\t%4.3f\t%4.3f\t%4.3f\t%4.3f\n",
                        input->worker[s]->label_string,
                        sd[s].locus->reference,
                        sd[s].locus->position,
                        sd[s].sample_points[p * 4],
                        sd[s].sample_points[p * 4 + 1],
                        sd[s].sample_points[p * 4 + 2],
                        sd[s].sample_points[p * 4 + 3]);
            }
            */

        }

        refresh_locus(input, s, &sd[s]);
    }

}




/* 
   input defines one range of loci for the worker to work on, across
   all samples and all pairings.  the 'worker' field of the input is a
   posterior_wrapper that holds its model parameters.

   [beg[s], end[s]) defines the workload.  it is assumed that, for any
   locus in each of these S ranges, and for any pairing of samples
   {s1, s2}, the union of loci in [beg[s1], end[s1]) and [beg[s2],
   end[s2]) define the set of locus pairs.  Any loci present in just
   one of the two samples will be filled in with the null locus for
   the other sample.

   cur[s] point to the current locus in [beg[s], end[s]) that is
   parsed and sampled.  All loci before cur[s] have already been
   processed.  The cur[s] loci are not necessarily the same locus, but
   the minimum among them defines the conceptual 'next' locus to
   generate distance comparisons for.

   In the event that [beg[s], end[s]) is an empty range for one or
   more 's', the null locus will be substituted for that sample.

   Once all pairwise distances have been generated for the 'next'
   locus, prints out compositions for all samples having data for this
   'next' locus.  Then, for that subset of samples, move cur[s], and
   refresh all dependent data structures.

   In order to progressively consume the loci in the range, there is a
   global_s, pointing to the least locus across all samples.
 */

void *dist_worker(void *args)
{
    dist_worker_input *input = static_cast<dist_worker_input *>(args);

    char *out_dist_buf = input->out_dist;
    char *out_comp_buf = input->out_comp;
    char *out_vcf_buf = input->out_vcf;

    if (out_dist_buf != NULL)
    {
        out_dist_buf[0] = '\0';
    }
    if (out_comp_buf != NULL)
    {
        out_comp_buf[0] = '\0';
    }
    if (out_vcf_buf != NULL)
    {
        out_vcf_buf[0] = '\0';
    }

    sample_details *samples = new sample_details[input->num_samples];
    size_t global_s = 0;

    for (size_t s = 0; s != input->num_samples; ++s)
    {
        samples[s].locus = NULL;
        samples[s].current = input->beg[s];
        samples[s].sample_points = new double[input->worker[s]->final_num_points * 4];
    }


    // before main loop, initialize loci and samples
    for (size_t s = 0; s != input->num_samples; ++s)
    {
        if (samples[s].current == input->end[s])
        {
            continue;
        }
        init_locus(input, s, &samples[s]);
    }

    global_s = init_global_and_next(input, samples);

    // main loop for computing pairwise distances
    // though slightly inefficient, just construct null_points here
    // !!! here, use worker[0] as a proxy.  this is a design flaw
    // owing to the fact that there are multiple workers used, all with the same
    // basic settings
    char null_algorithm[3];
    double *null_points = new double[input->worker[0]->final_num_points * 4];
    PileupSummary null_locus(0);

    input->worker[0]->params->pack(& null_locus.counts);
    input->worker[0]->sample(& null_locus, null_points, null_algorithm);

    while (samples[global_s].current != input->end[global_s])
    {
        // calculate pairwise distances with error bars for all pairs
        // for the current position.  Note that it is possible for
        // both samples in the pair to have no coverage, and the
        // distance metric will be the null distance distribution.
        // This is intentional.

        // print out distance quantiles for all sample pairs at the current locus
        out_dist_buf = next_distance_quantiles_aux(input, samples, global_s, 
                                                   null_points, out_dist_buf);

        // main loop for computing individual compositions
        // this is the last stop before we move forward
        advance_loci_aux(input, samples, & out_comp_buf, & out_vcf_buf);

        global_s = init_global_and_next(input, samples);
        
    }   

    for (size_t s = 0; s != input->num_samples; ++s)
    {
        delete samples[s].sample_points;
    }

    delete samples;
    delete null_points;

    pthread_exit((void*) 0);
}


