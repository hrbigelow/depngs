#include "dist_worker.h"
#include "pileup_tools.h"
#include "comp_functor.h"
#include "posterior.h"
// #include "integrands.h"
#include "error_estimate.h"

#include <cstdio>
#include <cstring>
#include <cassert>
#include <gsl/gsl_math.h>
#include <algorithm>

dist_worker_input::dist_worker_input(size_t num_samples,
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
                                     std::map<char const*, size_t, ltstr> * contig_order) :
    num_samples(num_samples), num_sample_pairs(num_sample_pairs),
    num_sample_point_pairings(num_sample_point_pairings),
    dist_quantiles(dist_quantiles), num_dist_quantiles(num_dist_quantiles),
    comp_quantiles(comp_quantiles), num_comp_quantiles(num_comp_quantiles),
    lattice(lattice),
    sampling_fallback_threshold(sampling_fallback_threshold),
    use_discrete(use_discrete), use_sampling(use_sampling),
    out_dist(out_dist), out_comp(out_comp), out_discomp(out_discomp),
    pair_sample1(pair_sample1),
    pair_sample2(pair_sample2), contig_order(contig_order)
{ 
    this->worker = new posterior_wrapper *[this->num_samples];
    this->beg = new std::vector<char *>::iterator[this->num_samples];
    this->end = new std::vector<char *>::iterator[this->num_samples];
    this->less_locus.contig_order = this->contig_order;
    this->equal_locus.contig_order = this->contig_order;
}


dist_worker_input::dist_worker_input() :
    worker(NULL), num_samples(0), num_sample_pairs(0),
    num_sample_point_pairings(0),
    dist_quantiles(NULL), num_dist_quantiles(0),
    comp_quantiles(NULL), num_comp_quantiles(0),
    lattice(NULL),
    use_discrete(false), use_sampling(false),
    beg(NULL),
    end(NULL),
    out_dist(NULL), out_comp(NULL),
    out_discomp(NULL),
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
}

// assumes the chromosome has a number
bool less_locus_position::operator()(char * locus_line1, char * locus_line2)
{
    char contig1[100];
    char contig2[100];
    size_t position1;
    size_t position2;
    sscanf(locus_line1, "%s\t%zu", contig1, &position1);
    sscanf(locus_line2, "%s\t%zu", contig2, &position2);

    if (strcmp(contig1, contig2) == 0) 
    { 
        return position1 < position2; 
    }
    else
    {
        std::map<char const*, size_t>::iterator it1, it2;
        it1 = this->contig_order->find(contig1);
        it2 = this->contig_order->find(contig2);
        assert(it1 != this->contig_order->end());
        assert(it2 != this->contig_order->end());
        return (*it1) < (*it2);
    }
}



bool equal_locus_position::operator()(char * locus_line1, char * locus_line2)
{
    char contig1[100];
    char contig2[100];
    size_t position1;
    size_t position2;
    sscanf(locus_line1, "%s\t%zu", contig1, &position1);
    sscanf(locus_line2, "%s\t%zu", contig2, &position2);
    return strcmp(contig1, contig2) == 0 && position1 == position2;
}


size_t distance_quantiles_locus_bytes(size_t num_quantiles)
{
    return 3 + 3 + 10 + 10 + (7 * num_quantiles) + (4 + num_quantiles);
}



// compute and print out distance quantiles, based on quasi-random
// pairings of two samplings for efficiency, the 'random' pairing is
// done simply by cycling through both sets of sample points, but
// starting in the middle for the second set.
// return the next write position
char * print_distance_quantiles(double const* points1,
                                double const* points2,
                                char const* contig,
                                size_t position,
                                dist_worker_input * wi,
                                char * out_dist_buf,
                                size_t pair_index)
{
    size_t s1 = wi->pair_sample1[pair_index];
    size_t s2 = wi->pair_sample2[pair_index];

    size_t num_random_pairs = wi->num_sample_point_pairings;
    size_t num_sample_points1 = wi->worker[s1]->final_num_points;
    size_t num_sample_points2 = wi->worker[s2]->final_num_points;
    double * square_dist = new double[num_random_pairs];

    for (size_t p = 0; p != num_random_pairs; ++p)
    {
        size_t pi1 = (p % num_sample_points1) * 4;
        size_t pi2 = ((p * 10) % num_sample_points2) * 4;
        square_dist[p] = 
            gsl_pow_2(points1[pi1] - points2[pi2])
            + gsl_pow_2(points1[pi1 + 1] - points2[pi2 + 1])
            + gsl_pow_2(points1[pi1 + 2] - points2[pi2 + 2])
            + gsl_pow_2(points1[pi1 + 3] - points2[pi2 + 3]);
    }

    // find the quantiles
    double * start = square_dist;
    double * end = square_dist + num_random_pairs;
    double * cut;
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


// compute and print out distance quantiles of the mutational distance
// at a locus between two samples.  The two samples have been
// evaluated at a set of discrete, pre-specified points, and the
// log(value) of the posterior is given for them.  the distance
// quantiles are calculated by first generating a set of weighted
// distances by taking the weighted cartesian product of all possible
// pairs that have a joint weight above a threshold.  The quantiles
// are then calculated from this set of weighted distances.
struct weighted_dist_pair
{
    double weight;
    double distance;
    weighted_dist_pair(double weight, double distance) : weight(weight), distance(distance) { }
    weighted_dist_pair() : weight(0), distance(0) { }
};

struct less_distance {
    bool operator()(weighted_dist_pair const& a, weighted_dist_pair const& b) const
    {
        return a.distance < b.distance;
    }
};



char * print_distance_quantiles_discrete(double const* posterior_values1,
                                         double const* posterior_values2,
                                         char const* contig,
                                         size_t position,
                                         dist_worker_input * wi,
                                         size_t pair_index,
                                         char * out_dist_buf)
{
    eval_dist_matrix * mat = wi->lattice;

    char const* sample_label1 = wi->worker[wi->pair_sample1[pair_index]]->label_string;
    char const* sample_label2 = wi->worker[wi->pair_sample2[pair_index]]->label_string;

    size_t num_filtered1, num_filtered2;
    size_t * index1 = new size_t[mat->num_points];
    size_t * index2 = new size_t[mat->num_points];

    // pack the indices
    size_t s1 = 0, s2 = 0;
    for (size_t p = 0; p != mat->num_points; ++p)
    {
        if (posterior_values1[p] >= mat->inclusion_threshold)
        {
            index1[s1] = p;
            ++s1;
        }
        if (posterior_values2[p] >= mat->inclusion_threshold)
        {
            index2[s2] = p;
            ++s2;
        }
    }

    num_filtered1 = s1;
    num_filtered2 = s2;
    size_t num_pairs = num_filtered1 * num_filtered2;
    assert(num_pairs > 0);

    // weighted_dist_pair * weighted_dist = new weighted_dist_pair[num_pairs];
    double * weights = new double[mat->num_distances];
    std::fill(weights, weights + mat->num_distances, 0.0);

    size_t p1, p2; // original discrete point indexes in [0, wi->num_eval_points)
    size_t i;
    double weight;
    double total_weight = 0.0;
    for (s1 = 0, i = 0; s1 != num_filtered1; ++s1, i += num_filtered2)
    {
        p1 = index1[s1];
        for (s2 = 0; s2 != num_filtered2; ++s2)
        {
            p2 = index2[s2];
            weight = posterior_values1[p1] * posterior_values2[p2];
            total_weight += weight;
            weights[mat->dist_index[p1 * mat->num_points + p2]] += weight;
        }
    }
    delete index1;
    delete index2;

    // normalize the weights, then calculate the cdf
    double cumul_weight = 0.0;
    for (size_t w = 0; w != mat->num_distances; ++w)
    {
        weights[w] /= total_weight;
        weights[w] += cumul_weight;
        cumul_weight = weights[w];
    }

    // do the specific exclusion test
    // in this test, we try to pre-emptively avoid the sort if we can determine that
    // the locus has zero values for all quantiles
    double highest_quantile = wi->dist_quantiles[wi->num_dist_quantiles - 1];
    // if (false)
    if (weights[0] > highest_quantile)
    {
        // do nothing
    }
    else
    {
        out_dist_buf += sprintf(out_dist_buf, "%s\t%s\t%s\t%Zu", 
                                sample_label1, sample_label2, contig, position);

        size_t w = 0;
        for (size_t q = 0; q != wi->num_dist_quantiles; ++q)
        {
            while (weights[w] < wi->dist_quantiles[q])
            {
                ++w;
            }
            out_dist_buf += sprintf(out_dist_buf, "\t%7.4f", mat->distances[w]);
        }    
        out_dist_buf += sprintf(out_dist_buf, "\n");
    }

    delete weights;

    return out_dist_buf;
}


struct sample_details
{
    PileupSummary * locus;
    bool is_next;
    double * sample_points;
    double * discrete_values;
    bool has_sample_points;
    std::vector<char *>::iterator current;
    char algorithm_used[10];
};

// print out distance quantiles for the next locus, for all pairs
char * next_distance_quantiles_aux(dist_worker_input * input, 
                                   sample_details * sd,
                                   size_t global_s,
                                   double const* null_points,
                                   double const* null_values,
                                   char * out_buf)
{
    if (out_buf == NULL)
    {
        return out_buf;
    }

    char contig[100];
    size_t position;
    sscanf(*sd[global_s].current, "%s\t%zu", contig, &position);

    for (size_t p = 0; p != input->num_sample_pairs; ++p)
    {
        size_t s1 = input->pair_sample1[p];
        size_t s2 = input->pair_sample2[p];
        sample_details & sd1 = sd[input->pair_sample1[p]];
        sample_details & sd2 = sd[input->pair_sample2[p]];

        bool do_sampling =
            (! input->use_discrete && input->use_sampling) // sampling is mandated
            || (input->use_sampling &&
                ((sd1.is_next && sd1.has_sample_points)
                 || (sd2.is_next && sd2.has_sample_points)
                 )
                );

        if (do_sampling)
        {
            // fill in points
            if (sd1.is_next && ! sd1.has_sample_points)
            {
                input->worker[s1]->sample(sd1.locus, sd1.sample_points, sd1.algorithm_used);
                sd1.has_sample_points = true;
            }
            if (sd2.is_next && ! sd2.has_sample_points)
            {
                input->worker[s2]->sample(sd2.locus, sd2.sample_points, sd2.algorithm_used);
                sd2.has_sample_points = true;
            }
            double const* points1 = sd1.is_next ? sd1.sample_points : null_points;
            double const* points2 = sd2.is_next ? sd2.sample_points : null_points;
            
            out_buf = 
                print_distance_quantiles(points1, points2, contig, position, input, out_buf, p);
        }
        else
        {
            // exclusively use discrete-based distance calculation
            double const* values1 = sd1.is_next ? sd1.discrete_values : null_values;
            double const* values2 = sd2.is_next ? sd2.discrete_values : null_values;
            
            out_buf = 
                print_distance_quantiles_discrete(values1, values2, contig, position, input, p, out_buf);
        }
        
    }
    return out_buf;
}


// initialize the locus defined by sd->current.
// assume sd->current is a valid iterator
void init_locus(dist_worker_input * input,
                size_t sample_id,
                sample_details * sd)
{
    size_t s = sample_id;
    sd->locus = new PileupSummary(0);
    sd->locus->load_line(*sd->current);
    sd->locus->parse(input->worker[s]->min_quality_score);
    input->worker[s]->posterior->ee->locus_data = & sd->locus->counts;
    input->worker[s]->params->pack(& sd->locus->counts);
                    
    // refresh sample points
    double max_posterior_val = 0.0;
    if (input->use_discrete)
    {
        input->worker[s]->values(input->lattice->points, input->lattice->num_points,
                                 sd->discrete_values);
        max_posterior_val = 
            * std::max_element(sd->discrete_values, sd->discrete_values + input->lattice->num_points);
    }
    bool do_sampling = 
        input->use_sampling
        && (! input->use_discrete
            || max_posterior_val < input->sampling_fallback_threshold);

    if (do_sampling)
    {
        input->worker[s]->sample(sd->locus, sd->sample_points, sd->algorithm_used);
        sd->has_sample_points = true;
    }
    else
    {
        sd->has_sample_points = false;
    }
}


// optionally destroy previous locus associated with sd,
// and create a new one, initializing it with sample points
// assume sd->locus and sd->current are initialized
void refresh_locus(dist_worker_input * input,
                   size_t sample_id,
                   sample_details * sd)
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
void advance_loci_aux(dist_worker_input * input,
                      sample_details * sd,
                      char ** out_comp_buf,
                      char ** out_discomp_buf)
{
    for (size_t s = 0; s != input->num_samples; ++s)
    {
        if (! sd[s].is_next)
        {
            continue;
        }

        if (*out_discomp_buf != NULL)
        {
            *out_discomp_buf =
                print_discrete_comp(sd[s].locus, 
                                    input->worker[s]->label_string,
                                    sd[s].discrete_values,
                                    input->lattice->num_points,
                                    input->lattice->max_points_to_print,
                                    input->lattice->min_value_to_print,
                                    *out_discomp_buf);
        }
        if (*out_comp_buf != NULL && sd[s].has_sample_points)
        {
            *out_comp_buf = 
                input->worker[s]->print_quantiles(sd[s].locus, sd[s].algorithm_used,
                                                  sd[s].sample_points, *out_comp_buf);
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

void * dist_worker(void * args)
{
    dist_worker_input * input = static_cast<dist_worker_input *>(args);

    char * out_dist_buf = input->out_dist;
    char * out_comp_buf = input->out_comp;
    char * out_discomp_buf = input->out_discomp;

    if (out_dist_buf != NULL)
    {
        out_dist_buf[0] = '\0';
    }
    if (out_comp_buf != NULL)
    {
        out_comp_buf[0] = '\0';
    }
    if (out_discomp_buf != NULL)
    {
        out_discomp_buf[0] = '\0';
    }

    sample_details * samples = new sample_details[input->num_samples];
    size_t global_s = 0;

    for (size_t s = 0; s != input->num_samples; ++s)
    {
        samples[s].locus = NULL;
        samples[s].current = input->beg[s];
        samples[s].sample_points = new double[input->worker[s]->final_num_points * 4];
        samples[s].discrete_values = new double[input->lattice->num_points];
    }


    // before main loop, initialize loci and samples
    for (size_t s = 0; s != input->num_samples; ++s)
    {
        if (samples[s].current == input->end[s])
        {
            continue;
        }
        init_locus(input, s, &samples[s]);

        // update global current to lesser position if needed
        // global_s should only be updated if both iterators are valid
        
        if (samples[s].current != input->end[s] &&
            samples[global_s].current != input->end[global_s] &&
            input->less_locus(*samples[s].current, *samples[global_s].current))
        {
            global_s = s;
        }
    }

    // initialize is_next
    for (size_t s = 0; s != input->num_samples; ++s)
    {
        samples[s].is_next = 
            (samples[s].current != input->end[s])
            && (samples[global_s].current != input->end[global_s])
            && input->equal_locus(*samples[s].current, *samples[global_s].current);
    }

    // main loop for computing pairwise distances
    // though slightly inefficient, just construct null_points here
    // !!! here, use worker[0] as a proxy.  this is a design flaw
    // owing to the fact that there are multiple workers used, all with the same
    // basic settings
    char null_algorithm[3];
    double * null_points = new double[input->worker[0]->final_num_points * 4];
    PileupSummary null_locus(0);

    input->worker[0]->params->pack(& null_locus.counts);
    input->worker[0]->sample(& null_locus, null_points, null_algorithm);

    double * null_values = new double[input->lattice->num_points];
    input->worker[0]->values(input->lattice->points, input->lattice->num_points, null_values);

    while (samples[global_s].current != input->end[global_s])
    {
        // calculate pairwise distances with error bars for all pairs
        // for the current position.  Note that it is possible for
        // both samples in the pair to have no coverage, and the
        // distance metric will be the null distance distribution.
        // This is intentional.

        // print out distance quantiles for all sample pairs at the current locus
        out_dist_buf = next_distance_quantiles_aux(input, samples, global_s, 
                                                   null_points, null_values, out_dist_buf);

        // main loop for computing individual compositions
        // this is the last stop before we move forward
        advance_loci_aux(input, samples, & out_comp_buf, & out_discomp_buf);

        for (size_t s = 0; s != input->num_samples; ++s)
        {
            if (samples[s].current != input->end[s]
                && samples[global_s].current != input->end[global_s]
                && input->less_locus(*samples[s].current, *samples[global_s].current))
            {
                global_s = s;
            }
        }
        
        // initialize is_next
        for (size_t s = 0; s != input->num_samples; ++s)
        {
            samples[s].is_next = 
                (samples[s].current != input->end[s]) 
                && (samples[global_s].current != input->end[global_s])
                && input->equal_locus(*samples[s].current, *samples[global_s].current);
        }

    }   

    for (size_t s = 0; s != input->num_samples; ++s)
    {
        delete samples[s].sample_points;
        delete samples[s].discrete_values;
    }

    delete samples;
    delete null_points;
    delete null_values;

    pthread_exit((void*) 0);
}
