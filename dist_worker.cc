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

dist_worker_input::dist_worker_input(size_t thread_num,
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
                                     char *out_dist,
                                     char *out_comp,
                                     char *out_vcf,
                                     char *out_indel_dist,
                                     size_t *pair_sample1,
                                     size_t *pair_sample2,
                                     std::map<const char *, size_t, ltstr> *contig_order) :
    thread_num(thread_num),
    num_samples(num_samples), num_sample_pairs(num_sample_pairs),
    num_sample_point_pairings(num_sample_point_pairings),
    dist_quantiles(dist_quantiles), num_dist_quantiles(num_dist_quantiles),
    comp_quantiles(comp_quantiles), num_comp_quantiles(num_comp_quantiles),
    min_high_conf_dist(min_high_conf_dist),
    prelim_num_points(prelim_num_points),
    prelim_quantile(prelim_quantile),
    final_num_points(final_num_points),
    out_dist(out_dist), out_comp(out_comp), out_vcf(out_vcf), out_indel_dist(out_indel_dist),
    pair_sample1(pair_sample1), pair_sample2(pair_sample2), 
    contig_order(contig_order)
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
    min_high_conf_dist(0),
    beg(NULL), end(NULL),
    out_dist(NULL), out_comp(NULL), out_vcf(NULL), out_indel_dist(NULL),
    pair_sample1(NULL), pair_sample2(NULL),
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




// compute distance quantiles by generating random pairs of points
// from the two sets of sample points
void compute_dist_quantiles(const double *points1,
                            const double *points2,
                            size_t num_dims,
                            size_t num_sample_points,
                            size_t num_random_pairs,
                            gsl_rng *randgen,
                            double *dist_quantiles,
                            size_t num_dist_quantiles,
                            double *dist_quantile_values)
{
    double 
        *square_dist = new double[num_random_pairs], 
        *sd = square_dist, 
        *sde = sd + num_random_pairs;

    const double *p1, *p2, *p1e;
    for (size_t p = 0; p != num_random_pairs; ++p, ++sd)
    {
        p1 = points1 + (p % num_sample_points) * num_dims, p1e = p1 + num_dims;
        p2 = points2 + gsl_rng_uniform_int(randgen, num_sample_points) * num_dims;
        *sd = 0;
        while (p1 != p1e) *sd += gsl_pow_2(*p1++ - *p2++);
    }

    // find the quantiles
    double *cut, qval;

    sd = square_dist;
    for (size_t q = 0; q != num_dist_quantiles; ++q)
    {
        cut = square_dist + static_cast<size_t>(std::round(dist_quantiles[q] * num_random_pairs));
        std::nth_element(sd, cut, sde);
        qval = (cut == sde) ? -1.0 : sqrt(*cut);
        sd = cut;
        dist_quantile_values[q] = qval;
    }
    delete[] square_dist;
}


struct indel_event
{
    unsigned count1, count2;
    const char *seq;
    bool is_insertion;
};


// compute and print out distance quantiles, based on quasi-random
// pairings of two samplings for efficiency, the 'random' pairing is
// done simply by cycling through both sets of sample points, but
// starting in the middle for the second set.
// return the next write position.
char *print_distance_quantiles(const char *contig,
                               size_t position,
                               dist_worker_input *wi,
                               size_t pair_index,
                               double *dist_quantile_values,
                               sample_details *sd,
                               char *out_dist_buf)
{
    size_t s1 = wi->pair_sample1[pair_index], s2 = wi->pair_sample2[pair_index];

    out_dist_buf += sprintf(out_dist_buf, "%s\t%s\t%s\t%Zu", 
                            wi->worker[s1]->label_string,
                            wi->worker[s2]->label_string,
                            contig, position);
    
    for (size_t q = 0; q != wi->num_dist_quantiles; ++q)
    {
        out_dist_buf += sprintf(out_dist_buf, "\t%07.4f", dist_quantile_values[q]);
    }

    if (sd)
        out_dist_buf += sprintf(out_dist_buf, "\t%Zu\t%s\t%s\t%Zu\t%s\t%s",
                                sd[s1].locus->read_depth,
                                sd[s1].locus->bases_raw,
                                sd[s1].locus->quality_codes,
                                sd[s2].locus->read_depth,
                                sd[s2].locus->bases_raw,
                                sd[s2].locus->quality_codes);

    out_dist_buf += sprintf(out_dist_buf, "\n");

    return out_dist_buf;
}



char *print_indel_distance_quantiles(const char *contig,
                                     size_t position,
                                     dist_worker_input *wi,
                                     size_t pair_index,
                                     double *dist_quantile_values,
                                     indel_event *events,
                                     size_t num_events,
                                     sample_details *sd,
                                     char *out_dist_buf)
{
    size_t s1 = wi->pair_sample1[pair_index], s2 = wi->pair_sample2[pair_index];

    out_dist_buf += sprintf(out_dist_buf, "%s\t%s\t%s\t%Zu", 
                            wi->worker[s1]->label_string,
                            wi->worker[s2]->label_string,
                            contig, position);
    
    for (size_t q = 0; q != wi->num_dist_quantiles; ++q)
    {
        out_dist_buf += sprintf(out_dist_buf, "\t%.4f", dist_quantile_values[q]);
    }

    // now print the indel event summary
    indel_event *eb = events, *ee = eb + num_events;
    while (eb != ee) 
    {
        out_dist_buf += sprintf(out_dist_buf, "%c%i", (eb == events ? '\t' : ','), eb->count1);
        eb++;
    }
    eb = events;
    while (eb != ee) 
    {
        out_dist_buf += sprintf(out_dist_buf, "%c%i", (eb == events ? '\t' : ','), eb->count2);
        eb++;
    }
    eb = events;
    while (eb != ee) 
    {
        out_dist_buf += sprintf(out_dist_buf, "%c%c%s", 
                                (eb == events ? '\t' : ','), 
                                (eb->seq ? (eb->is_insertion ? '+' : '-') : '@'),
                                (eb->seq ? eb->seq : ""));
        eb++;
    }

    if (sd)
        out_dist_buf += sprintf(out_dist_buf, "\t%Zu\t%s\t%s\t%Zu\t%s\t%s",
                                sd[s1].locus->read_depth,
                                sd[s1].locus->bases_raw,
                                sd[s1].locus->quality_codes,
                                sd[s2].locus->read_depth,
                                sd[s2].locus->bases_raw,
                                sd[s2].locus->quality_codes);

    out_dist_buf += sprintf(out_dist_buf, "\n");

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
    const double *mode1, *mode2, *end1;
    double *dist_quantile_values = new double[input->num_dist_quantiles];
    double min_dist_squared = input->min_high_conf_dist * input->min_high_conf_dist;

    sscanf(*sd[global_s].current, "%s\t%zu", contig, &position);

    for (size_t p = 0; p != input->num_sample_pairs; ++p)
    {
        size_t s1 = input->pair_sample1[p], s2 = input->pair_sample2[p];
        sample_details *sd1 = &sd[s1], *sd2 = &sd[s2];

        if (input->min_high_conf_dist > 0 && ! (sd1->is_next && sd2->is_next))
            continue;
        
        if (sd1->is_next && ! sd1->mode_computed) (input->worker[s1]->find_mode(), sd1->mode_computed = true);
        if (sd2->is_next && ! sd2->mode_computed) (input->worker[s2]->find_mode(), sd2->mode_computed = true);

        mode1 = sd1->is_next ? input->worker[s1]->mode_point : NULL_MODE, end1 = mode1 + 4;
        mode2 = sd2->is_next ? input->worker[s2]->mode_point : NULL_MODE;

        float mode_square_dist = 0;
        while (mode1 != end1)
            mode_square_dist += gsl_pow_2(*mode1++ - *mode2++);

        if (mode_square_dist < min_dist_squared)
        {
            // modes don't differ enough
            continue;
        }

        // preliminary sampling
        if (sd1->is_next && ! sd1->num_sample_points)
        {
            input->worker[s1]->tune(sd1);
            input->worker[s1]->sample(sd1, input->prelim_num_points);
            sd1->num_sample_points = input->prelim_num_points;
        }
        if (sd2->is_next && ! sd2->num_sample_points)
        {
            input->worker[s2]->tune(sd2);
            input->worker[s2]->sample(sd2, input->prelim_num_points);
            sd2->num_sample_points = input->prelim_num_points;
        }

        compute_dist_quantiles(sd1->is_next ? sd1->sample_points : null_points,
                               sd2->is_next ? sd2->sample_points : null_points,
                               4,
                               input->prelim_num_points,
                               input->prelim_num_points * 10, // ad-hoc
                               input->randgen,
                               &input->prelim_quantile, 1,
                               dist_quantile_values);
        
        if (dist_quantile_values[0] < input->min_high_conf_dist)
        {
            // after prelim testing, still not enough difference
            continue;
        }
        
        if (sd1->num_sample_points == input->prelim_num_points)
        {
            input->worker[s1]->sample(sd1, input->final_num_points);
            sd1->num_sample_points = input->final_num_points;
        }
        if (sd2->num_sample_points == input->prelim_num_points)
        {
            input->worker[s2]->sample(sd2, input->final_num_points);
            sd2->num_sample_points = input->final_num_points;
        }

        compute_dist_quantiles(sd1->is_next ? sd1->sample_points : null_points,
                               sd2->is_next ? sd2->sample_points : null_points,
                               4,
                               input->final_num_points,
                               input->num_sample_point_pairings,
                               input->randgen,
                               input->dist_quantiles,
                               input->num_dist_quantiles,
                               dist_quantile_values);
        
        if (dist_quantile_values[0] >= input->min_high_conf_dist)
        {
            out_buf = 
                print_distance_quantiles(contig, position, input, p, dist_quantile_values, NULL, out_buf);
        }
    }
    delete[] dist_quantile_values;
    return out_buf;
}

indel_event *set_indel_events_aux(CHAR_MAP::iterator id1,
                                  CHAR_MAP::iterator e1,
                                  CHAR_MAP::iterator id2,
                                  CHAR_MAP::iterator e2,
                                  bool is_insertion,
                                  indel_event *e)
{
    while (id1 != e1 || id2 != e2)
    {
        e->count1 = id1 != e1 && (id2 == e2 || strcmp((*id1).first, (*id2).first) <= 0) ? (*id1).second : 0;
        e->count2 = id2 != e2 && (id1 == e1 || strcmp((*id2).first, (*id1).first) <= 0) ? (*id2).second : 0;
        if (e->count1 != 0) { e->seq = (*id1).first; ++id1; }
        if (e->count2 != 0) { e->seq = (*id2).first; ++id2; }
        e->is_insertion = is_insertion;
        ++e;
    }
    return e;
}

// generate a tally of counts of each type of indel and return a
// allocated array of the counts
indel_event *count_indel_types(sample_details *sd1,
                            sample_details *sd2, 
                            size_t *num_counts)
{

    if (! (sd1->is_next && sd2->is_next))
    {
        *num_counts = 0;
        return NULL;
    }
    else
    {
        if (sd1->locus->position == 761957)
        {
            int i = 0;
            ++i;
        }

        *num_counts = 0;

        size_t max_possible = 
            sd1->locus->deletions.size() + sd1->locus->insertions.size()
            + sd2->locus->deletions.size() + sd2->locus->insertions.size()
            + 1;
        
        indel_event *events = new indel_event[max_possible], *e = events, *ee = events;

        CHAR_MAP::iterator id1, id2, e1, e2;

        id1 = sd1->locus->deletions.begin(), e1 = sd1->locus->deletions.end();
        id2 = sd2->locus->deletions.begin(), e2 = sd2->locus->deletions.end();
        e = set_indel_events_aux(id1, e1, id2, e2, false, e);
        
        id1 = sd1->locus->insertions.begin(), e1 = sd1->locus->insertions.end();
        id2 = sd2->locus->insertions.begin(), e2 = sd2->locus->insertions.end();
        e = set_indel_events_aux(id1, e1, id2, e2, true, e);

        // count numbers of indels
        unsigned nindel1 = 0, nindel2 = 0;
        while (ee != e) nindel1 += ee->count1, nindel2 += ee++->count2;

        // now count the non-indel 'events', by definition the
        // remaining reads that do not have any indels.  only count
        // this as an event if it occurs in at least one of the two
        // pairs.
        if ((e->count1 = sd1->locus->read_depth - nindel1) +
            (e->count2 = sd2->locus->read_depth - nindel2) > 0)
        {
            e->seq = NULL;
            ++e;
        }

        *num_counts = e - events;

        return events;
    }    
}






// print out all next distance quantiles for indels
char *next_indel_distance_quantiles_aux(dist_worker_input *input, 
                                        sample_details *sd,
                                        size_t global_s,
                                        char *out_buf)
{
    /*
      1.  count the indel types
      2.  allocate two buffers
      3.  populate each buffer with dirichlets
      4.  compute the dist quantiles
      5.  deallocate the buffers
      5.  print out suitably filtered distances
    */    

    if (out_buf == NULL) return out_buf;

    char contig[100];
    size_t position, num_events;
    double *dist_quantile_values = new double[input->num_dist_quantiles];
    
    sscanf(*sd[global_s].current, "%s\t%zu", contig, &position);

    for (size_t p = 0; p != input->num_sample_pairs; ++p)
    {
        size_t s1 = input->pair_sample1[p], s2 = input->pair_sample2[p];
        sample_details *sd1 = &sd[s1], *sd2 = &sd[s2];

        if (! (sd1->is_next && sd2->is_next) && input->min_high_conf_dist > 0) continue;
        
        indel_event *all_events = count_indel_types(sd1, sd2, &num_events);

        if (num_events >= 2) 
        {
            // need at least some indels, otherwise these loci don't differ
            // there will always be at least one event, the reads themselves.  (though the count may be zero)
            double *alpha1 = new double[num_events], *alpha2 = new double[num_events];
            for (size_t c = 0; c != num_events; ++c) 
            {
                alpha1[c] = all_events[c].count1 + 1;
                alpha2[c] = all_events[c].count2 + 1;
            }
        
            size_t bufsize = num_events * input->final_num_points;
            double *points1 = new double[bufsize], *p1 = points1, *pe1 = points1 + bufsize;
            double *points2 = new double[bufsize], *p2 = points2;

            while (p1 != pe1) 
            {
                gsl_ran_dirichlet(input->randgen, num_events, alpha1, p1);
                gsl_ran_dirichlet(input->randgen, num_events, alpha2, p2);
                p1 += num_events;
                p2 += num_events;
            }

            // compute distance quantiles
            compute_dist_quantiles(points1, points2, num_events, 
                                   input->final_num_points,
                                   input->num_sample_point_pairings,
                                   input->randgen,
                                   input->dist_quantiles,
                                   input->num_dist_quantiles,
                                   dist_quantile_values);

            // print out suitably filtered distances
            if (dist_quantile_values[0] >= input->min_high_conf_dist)
            {
                out_buf = 
                    print_indel_distance_quantiles(contig, position, input, p, 
                                                   dist_quantile_values, 
                                                   all_events, num_events, sd, out_buf);
            }

            delete[] points1;
            delete[] points2;
            delete[] alpha1;
            delete[] alpha2;
        }
        delete[] all_events;
    }
    delete[] dist_quantile_values;
    return out_buf;

}


/*
void print_indel_comparisons(dist_worker_input *wi,
                             sample_details *sd,
                             size_t global_s)
{

    std::map<unsigned, unsigned>::iterator dit1, e1, dit2, e2;
    CHAR_MAP::iterator iit;

    char contig[100];
    size_t position;
    sscanf(*sd[global_s].current, "%s\t%zu", contig, &position);

    for (size_t p = 0; p != wi->num_sample_pairs; ++p)
    {
        size_t s1 = wi->pair_sample1[p], s2 = wi->pair_sample2[p];
        sample_details *sd1 = &sd[s1], *sd2 = &sd[s2];

        unsigned max1 = (! sd1->is_next) || sd1->locus->deletions.empty() 
            ? 0 : (*sd1->locus->deletions.rbegin()).first;

        unsigned max2 = (! sd2->is_next) || sd2->locus->deletions.empty()
            ? 0 : (*sd2->locus->deletions.rbegin()).first;

        unsigned max = max1 > max2 ? max1 : max2;

        if (max != 0)
        {
        
            printf("%s\t%s\t%s\t%Zu\t%Zu\t%Zu",
                   wi->worker[s1]->label_string,
                   wi->worker[s2]->label_string,
                   contig,
                   position,
                   sd1->locus->read_depth,
                   sd2->locus->read_depth);



            dit1 = sd1->locus->deletions.begin(), e1 = sd1->locus->deletions.end();
            dit2 = sd2->locus->deletions.begin(), e2 = sd2->locus->deletions.end();

            unsigned *cnt1 = new unsigned[max * 10], *cnt2 = new unsigned[max * 10];

            unsigned del = 0;
            while (del <= max)
            {
                cnt1[del] = (dit1 != e1 && (*dit1).first == del ? (*dit1++).second : 0);
                cnt2[del] = (dit2 != e2 && (*dit2).first == del ? (*dit2++).second : 0);
                ++del;
            }

            for (size_t d = 0; d != del; ++d)
                printf("%c%u", (d == 0 ? '\t' : ','), cnt1[d]);

            for (size_t d = 0; d != del; ++d)
                printf("%c%u", (d == 0 ? '\t' : ','), cnt2[d]);

            printf("\n");

            delete[] cnt1;
            delete[] cnt2;
        }
    }
}
*/


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
    sd->locus = new PileupSummary;
    sd->locus->load_line(*sd->current);
    sd->locus->parse(input->worker[s]->min_quality_score);
    sd->num_sample_points = 0;

    input->worker[s]->model->locus_data = &sd->locus->counts;
    input->worker[s]->params->pack(&sd->locus->counts);
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
        if (! sd[s].is_next) continue;

        if (*out_comp_buf != NULL && sd[s].num_sample_points)
        {
            *out_comp_buf = input->worker[s]->print_quantiles(&sd[s], *out_comp_buf);

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
    char *out_indel_dist_buf = input->out_indel_dist;

    if (out_dist_buf != NULL) out_dist_buf[0] = '\0';
    if (out_comp_buf != NULL) out_comp_buf[0] = '\0';
    if (out_vcf_buf != NULL) out_vcf_buf[0] = '\0';
    if (out_indel_dist_buf != NULL) out_indel_dist_buf[0] = '\0';

    sample_details *samples = new sample_details[input->num_samples];
    size_t global_s = 0;

    for (size_t s = 0; s != input->num_samples; ++s)
    {
        samples[s].locus = NULL;
        samples[s].current = input->beg[s];
        samples[s].sample_points = new double[input->final_num_points * 4];
        samples[s].num_sample_points = 0;
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
    PileupSummary null_locus;
    input->worker[0]->model->locus_data = &null_locus.counts;
    input->worker[0]->params->pack(& null_locus.counts);
    sample_details null_sd;
    null_sd.locus = &null_locus;
    null_sd.sample_points = new double[input->final_num_points * 4];
    null_sd.mode_computed = true;
    input->worker[0]->find_mode();
    input->worker[0]->tune(&null_sd);
    input->worker[0]->sample(&null_sd, input->final_num_points);

    while (samples[global_s].current != input->end[global_s])
    {
        // print out distance quantiles for all sample pairs at the current locus
        out_dist_buf = next_distance_quantiles_aux(input, samples, global_s, 
                                                   null_sd.sample_points,
                                                   out_dist_buf);

        out_indel_dist_buf = 
            next_indel_distance_quantiles_aux(input, samples, global_s, out_indel_dist_buf);

        // main loop for computing individual compositions
        // this is the last stop before we move forward
        advance_loci_aux(input, samples, & out_comp_buf, & out_vcf_buf);

        // print_indel_comparisons(input, samples, global_s);

        global_s = init_global_and_next(input, samples);
        
    }   

    delete[] null_sd.sample_points;
    for (size_t s = 0; s != input->num_samples; ++s)
    {
        delete[] samples[s].sample_points;
    }

    delete[] samples;

    pthread_exit((void*) 0);
}


