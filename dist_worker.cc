#include "dist_worker.h"
#include "pileup_tools.h"
#include "comp_worker.h"
#include "defs.h"
#include "metropolis_sampling.h"

#include <gsl/gsl_math.h>
#include <algorithm>

extern "C" {
#include "locus.h"
}

#define MIN(a, b) ((a) < (b) ? (a) : (b))


// compute distance quantiles by generating random pairs of points
// from the two sets of sample points
void compute_dist_quantiles(const double *points1,
                            const double *points2,
                            size_t n_dims,
                            size_t n_sample_points,
                            double *square_dist_buf,
                            size_t n_random_pairs,
                            gsl_rng *randgen,
                            const double *dist_quantiles,
                            size_t n_dist_quantiles,
                            double *dist_quantile_values)
{
    double 
        *sd = square_dist_buf, 
        *sde = sd + n_random_pairs;

    const double *p1, *p2, *p1e;
    for (size_t p = 0; p != n_random_pairs; ++p, ++sd)
    {
        p1 = points1 + (p % n_sample_points) * n_dims, p1e = p1 + n_dims;
        p2 = points2 + gsl_rng_uniform_int(randgen, n_sample_points) * n_dims;
        *sd = 0;
        while (p1 != p1e) *sd += gsl_pow_2(*p1++ - *p2++);
    }

    // find the quantiles
    double *cut, qval;

    sd = square_dist_buf;
    for (size_t q = 0; q != n_dist_quantiles; ++q)
    {
        cut = square_dist_buf + (size_t)ceil(dist_quantiles[q] * n_random_pairs);
        std::nth_element(sd, cut, sde);
        qval = (cut == sde) ? -1.0 : sqrt(*cut);
        sd = cut;
        dist_quantile_values[q] = qval;
    }
}


// print out distance quantiles, based on quasi-random pairings of two
// samplings for efficiency, the 'random' pairing is done simply by
// cycling through both sets of sample points, but starting in the
// middle for the second set.  return the next write position.
char *print_distance_quantiles(const char *contig,
                               size_t position,
                               dist_worker_input *wi,
                               size_t pair_index,
                               double *dist_quantile_values,
                               locus_sampling *sd,
                               char *out_dist_buf)
{
    size_t s1 = wi->pair_sample1[pair_index], s2 = wi->pair_sample2[pair_index];

    out_dist_buf += sprintf(out_dist_buf, "%s\t%s\t%s\t%Zu", 
                            wi->sample_atts[s1].label_string,
                            wi->sample_atts[s2].label_string,
                            contig, position);
    
    for (size_t q = 0; q != wi->pset->n_dist_quantiles; ++q)
       out_dist_buf += sprintf(out_dist_buf, "\t%7.4f", dist_quantile_values[q]);

    if (wi->print_pileup_fields)
        out_dist_buf += sprintf(out_dist_buf, "\t%Zu\t%s\t%s\t%Zu\t%s\t%s",
                                sd[s1].locus.read_depth,
                                sd[s1].locus.bases_raw.buf,
                                sd[s1].locus.quality_codes.buf,
                                sd[s2].locus.read_depth,
                                sd[s2].locus.bases_raw.buf,
                                sd[s2].locus.quality_codes.buf);

    sd[s1].dist_printed = 1;
    sd[s2].dist_printed = 1;

    out_dist_buf += sprintf(out_dist_buf, "\n");

    return out_dist_buf;
}

/* summarizes the counts of each indel in the pair of samples
   together, occurring at a particular locus.  the indel event is
   associated with a single base locus, even though, for example, a
   deletion may span multiple loci.  by convention, the locus that
   occurs just before the inserted or deleted dna is the locus
   associated with the indel event.  */
struct indel_event
{
    unsigned count1, count2;
    const char *seq; // the sequence that is either deleted or inserted
    bool is_insertion; // whether or not this is an insertion
};


char *print_indel_distance_quantiles(const char *contig,
                                     size_t position,
                                     dist_worker_input *wi,
                                     size_t pair_index,
                                     double *dist_quantile_values,
                                     indel_event *events,
                                     size_t n_events,
                                     locus_sampling *sd,
                                     char *out_dist_buf)
{
    size_t s1 = wi->pair_sample1[pair_index], s2 = wi->pair_sample2[pair_index];

    out_dist_buf += sprintf(out_dist_buf, "%s\t%s\t%s\t%Zu", 
                            wi->sample_atts[s1].label_string,
                            wi->sample_atts[s2].label_string,
                            contig, position);
    
    for (size_t q = 0; q != wi->pset->n_dist_quantiles; ++q)
    {
        out_dist_buf += sprintf(out_dist_buf, "\t%.4f", dist_quantile_values[q]);
    }

    // now print the indel event summary
    indel_event *eb = events, *ee = eb + n_events;
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
                                sd[s1].locus.read_depth,
                                sd[s1].locus.bases_raw.buf,
                                sd[s1].locus.quality_codes.buf,
                                sd[s2].locus.read_depth,
                                sd[s2].locus.bases_raw.buf,
                                sd[s2].locus.quality_codes.buf);

    out_dist_buf += sprintf(out_dist_buf, "\n");

    return out_dist_buf;

}


/* print out distance quantiles for the next locus, for all pairs.
   also generates sample points for each sample as needed, both for
   the preliminary test and more points for the final test */
char *next_distance_quantiles_aux(dist_worker_input *input, 
                                  locus_sampling *sd,
                                  size_t gs,
                                  const double *null_points,
                                  char *out_buf)
{
    if (! out_buf) return out_buf;

    double estimated_mean[2][NUM_NUCS], proposal_alpha[2][NUM_NUCS];
    locus_sampling *sdpair[2];
    size_t cumul_aoff[2];
    size_t p, i;

    for (p = 0; p != input->n_sample_pairs; ++p)
    {
        sdpair[0] = &sd[input->pair_sample1[p]];
        sdpair[1] = &sd[input->pair_sample2[p]];

        if (input->min_high_conf_dist > 0 && ! (sdpair[0]->is_next && sdpair[1]->is_next))
            continue;
        
        /* preliminary sampling */
        for (i = 0; i != 2; ++i)
        {
            if (sdpair[i]->is_next && ! sdpair[i]->n_sample_points)
            {
                cumul_aoff[i] = 
                    tune_proposal(&sdpair[i]->locus.counts, 
                                  input->pset,
                                  proposal_alpha[i], estimated_mean[i], 
                                  sdpair[i]->sample_points);
                
                metropolis_sampling(0, input->prelim_n_points, 
                                    &sdpair[i]->locus.counts,
                                    input->pset->logu, proposal_alpha[i], 
                                    input->pset->prior_alpha,
                                    cumul_aoff[i], 
                                    sdpair[i]->sample_points);
                
                // input->worker[s1]->tune(sd1, estimated_mean1);
                // input->worker[s1]->sample(sd1, estimated_mean1, input->prelim_n_points);
                sdpair[i]->n_sample_points = input->prelim_n_points;
            }
        }
            
        compute_dist_quantiles(sdpair[0]->is_next ? sdpair[0]->sample_points : null_points,
                               sdpair[1]->is_next ? sdpair[1]->sample_points : null_points,
                               4,
                               input->prelim_n_points, input->square_dist_buf,
                               MIN(input->n_sample_point_pairings,
                                   input->prelim_n_points * 10), // ad-hoc
                               input->randgen,
                               &input->prelim_quantile, 1, input->dist_quantile_values);
        
        /* after prelim testing, still not enough difference */
        if (input->dist_quantile_values[0] < input->min_high_conf_dist)
            continue;
        
        for (i = 0; i != 2; ++i)
        {
            if (sdpair[i]->n_sample_points == input->prelim_n_points)
            {
                metropolis_sampling(sdpair[i]->n_sample_points, 
                                    input->pset->final_n_points, 
                                    &sdpair[i]->locus.counts,
                                    input->pset->logu, proposal_alpha[i], 
                                    input->pset->prior_alpha,
                                    cumul_aoff[i], 
                                    sdpair[i]->sample_points);
                
                // input->worker[s1]->sample(sdpair[i], estimated_mean1, input->final_n_points);
                sdpair[i]->n_sample_points = input->pset->final_n_points;
            }
        }

        compute_dist_quantiles(sdpair[0]->is_next ? sdpair[0]->sample_points : null_points,
                               sdpair[1]->is_next ? sdpair[1]->sample_points : null_points,
                               4,
                               input->pset->final_n_points,
                               input->square_dist_buf,
                               input->n_sample_point_pairings,
                               input->randgen,
                               input->pset->dist_quantiles,
                               input->pset->n_dist_quantiles,
                               input->dist_quantile_values);
        
        if (input->dist_quantile_values[0] >= input->min_high_conf_dist)
            out_buf = 
                print_distance_quantiles(sd[gs].locus.reference, 
                                         sd[gs].locus.position, input, 
                                         p, input->dist_quantile_values, 
                                         sd, out_buf);
    }
    return out_buf;
}


/*  id1, e1 is the range over the first sample's insertions (or
    deletions), id2, e2 is the range over the second sample's
    insertions (or deletions). this function is called once for
    insertions, once for deletions, on each locus.  initializes as
    many indel_event's as needed.  automatically detects co-occurring
    insertions (deletions) and singly-occuring ones.  */
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

/* generate a tally of counts of each type of indel and return a
   allocated array of the counts */
indel_event *count_indel_types(locus_sampling *sd1,
                               locus_sampling *sd2, 
                               size_t *n_counts)
{

    if (! (sd1->is_next && sd2->is_next))
    {
        *n_counts = 0;
        return NULL;
    }
    else
    {
        // if (sd1->locus.position == 761957)
        // {
        //     int i = 0;
        //     ++i;
        // }

        *n_counts = 0;

        size_t max_possible = 
            sd1->locus.deletions.size() + sd1->locus.insertions.size()
            + sd2->locus.deletions.size() + sd2->locus.insertions.size()
            + 1;
        
        indel_event *events = new indel_event[max_possible], *e = events, *ee = events;

        CHAR_MAP::iterator id1, id2, e1, e2;

        id1 = sd1->locus.deletions.begin(), e1 = sd1->locus.deletions.end();
        id2 = sd2->locus.deletions.begin(), e2 = sd2->locus.deletions.end();
        e = set_indel_events_aux(id1, e1, id2, e2, false, e);
        
        id1 = sd1->locus.insertions.begin(), e1 = sd1->locus.insertions.end();
        id2 = sd2->locus.insertions.begin(), e2 = sd2->locus.insertions.end();
        e = set_indel_events_aux(id1, e1, id2, e2, true, e);

        // count numbers of indels
        unsigned nindel1 = 0, nindel2 = 0;
        while (ee != e) nindel1 += ee->count1, nindel2 += ee++->count2;

        // now count the non-indel 'events', by definition the
        // remaining reads that do not have any indels.  only count
        // this as an event if it occurs in at least one of the two
        // pairs.
        if ((e->count1 = sd1->locus.read_depth - nindel1) +
            (e->count2 = sd2->locus.read_depth - nindel2) > 0)
        {
            e->seq = NULL;
            ++e;
        }

        *n_counts = e - events;

        return events;
    }    
}


// print out all next distance quantiles for indels
char *next_indel_distance_quantiles_aux(dist_worker_input *input, 
                                        locus_sampling *sd,
                                        size_t gs,
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
    size_t position, n_events;
    
    sscanf(sd[gs].current, "%s\t%zu", contig, &position);

    for (size_t p = 0; p != input->n_sample_pairs; ++p)
    {
        size_t s1 = input->pair_sample1[p], s2 = input->pair_sample2[p];
        locus_sampling *sd1 = &sd[s1], *sd2 = &sd[s2];

        if (! (sd1->is_next && sd2->is_next) && input->min_high_conf_dist > 0) continue;
        
        indel_event *all_events = count_indel_types(sd1, sd2, &n_events);

        if (n_events >= 2) 
        {
            // need at least some indels, otherwise these loci don't differ
            // there will always be at least one event, the reads themselves.  (though the count may be zero)
            double *alpha1 = new double[n_events], *alpha2 = new double[n_events];
            for (size_t c = 0; c != n_events; ++c) 
            {
                alpha1[c] = all_events[c].count1 + 1;
                alpha2[c] = all_events[c].count2 + 1;
            }
        
            size_t bufsize = n_events * input->pset->final_n_points;
            double *points1 = new double[bufsize], *p1 = points1, *pe1 = points1 + bufsize;
            double *points2 = new double[bufsize], *p2 = points2;

            while (p1 != pe1) 
            {
                gsl_ran_dirichlet(input->randgen, n_events, alpha1, p1);
                gsl_ran_dirichlet(input->randgen, n_events, alpha2, p2);
                p1 += n_events;
                p2 += n_events;
            }

            // compute distance quantiles
            compute_dist_quantiles(points1,
                                   points2, 
                                   n_events, 
                                   input->pset->final_n_points,
                                   input->square_dist_buf,
                                   input->n_sample_point_pairings,
                                   input->randgen,
                                   input->pset->dist_quantiles,
                                   input->pset->n_dist_quantiles,
                                   input->dist_quantile_values);

            // print out suitably filtered distances
            if (input->dist_quantile_values[0] >= input->min_high_conf_dist)
            {
                out_buf = 
                    print_indel_distance_quantiles(contig, position, input, p, 
                                                   input->dist_quantile_values, 
                                                   all_events, n_events, sd, out_buf);
            }

            delete[] points1;
            delete[] points2;
            delete[] alpha1;
            delete[] alpha2;
        }
        delete[] all_events;
    }
    return out_buf;

}


// find gs, and init the is_next field for all samples
// this is run
size_t init_global_and_next(dist_worker_input *input, locus_sampling *samples)
{
    size_t gs = 0;
    for (size_t s = 0; s != input->n_samples; ++s)
    {
        if (samples[s].current != samples[s].end
            && (samples[gs].current == samples[gs].end
                || cmp_pair_ordering(&samples[s].locus_ord, 
                                     &samples[gs].locus_ord) < 0)
            )
            gs = s;
    }
    bool global_at_end = (samples[gs].current == samples[gs].end);
    for (size_t s = 0; s != input->n_samples; ++s)
        samples[s].is_next = 
            samples[s].current != samples[s].end
            && (! global_at_end)
            && cmp_pair_ordering(&samples[s].locus_ord,
                                 &samples[gs].locus_ord) == 0;

    return gs;
}


/* receives a certain number of in_bufs and a certain number of
   out_bufs.  par (cast to struct dist_worker_input) tells dist_worker
   how many input and output buffers to expect, and how to use them.
   
   there is one struct locus_sampling for each input.  it's current
   field points to the current line being processed. 'gs' is a single
   index indicating the sample with the lowest 'current' among all of
   them.  it is this position that must be fully processed before any
   samples may advance.

   any sample missing the locus defined by sample[gs].current has
   null_sd substituted for it.
 */
void dist_worker(void *par, const struct managed_buf *in_bufs,
                 struct managed_buf *out_bufs)
{
    struct dist_worker_input *dw = (struct dist_worker_input *)par;

    unsigned i = 0;
    struct managed_buf 
        *dist_buf = dw->do_dist ? &out_bufs[i++] : NULL,
        *comp_buf = dw->do_comp ? &out_bufs[i++] : NULL,
        *indel_buf = dw->do_indel ? &out_bufs[i++] : NULL;

    char
        *dist_ptr = dist_buf ? dist_buf->buf : NULL,
        *comp_ptr = comp_buf ? comp_buf->buf : NULL,
        *indel_ptr = indel_buf ? indel_buf->buf : NULL;
    
    if (dist_ptr) dist_ptr[0] = '\0';
    if (comp_ptr) comp_ptr[0] = '\0';
    if (indel_ptr) indel_ptr[0] = '\0';
    
    /* struct locus_sampling *samples = (struct locus_sampling *)
       malloc(dw->n_samples * sizeof(struct locus_sampling)); */
    /* 'new' is needed here because the struct contains a
       PileupSummary, which has a std::map. */
    struct locus_sampling *samples = new struct locus_sampling[dw->n_samples];

    size_t gs = 0, s;

    for (s = 0; s != dw->n_samples; ++s)
    {    
        samples[s].locus = PileupSummary();
        samples[s].is_next = false;
        samples[s].dist_printed = 0;
        samples[s].sample_points = 
            (double *)malloc(dw->pset->final_n_points * 4 * sizeof(double));
        samples[s].n_sample_points = 0;
        samples[s].autocor_offset = 0;
        samples[s].current = in_bufs[s].buf;
        samples[s].end = in_bufs[s].buf + in_bufs[s].size;
        samples[s].locus_ord = min_pair_ord;
    }

    // before main loop, initialize loci and samples
    for (s = 0; s != dw->n_samples; ++s)
    {
        if (samples[s].current == samples[s].end) continue;
        init_pileup_locus(&dw->sample_atts[s].nuc_stats, 
                          dw->pset->min_quality_score, &samples[s]);
    }

    gs = init_global_and_next(dw, samples);

    // main loop for computing pairwise distances
    // though slightly inefficient, just construct null_points here
    // !!! here, use worker[0] as a proxy.  this is a design flaw
    // owing to the fact that there are multiple workers used, all with the same
    // basic settings
    locus_sampling null_sd;
    null_sd.locus = PileupSummary();
    null_sd.sample_points = 
        (double *)malloc(dw->pset->final_n_points * 4 * sizeof(double));
    nucleotide_stats_pack(&dw->sample_atts[0].nuc_stats, &null_sd.locus.counts);
    
    double estimated_mean[NUM_NUCS], proposal_alpha[NUM_NUCS];;
    size_t cumul_aoff = tune_proposal(&null_sd.locus.counts,
                                      dw->pset, 
                                      proposal_alpha, estimated_mean,
                                      null_sd.sample_points);
    
    // dw->worker[0]->tune(&null_sd, estimated_mean);
    metropolis_sampling(0, dw->pset->final_n_points,
                        &null_sd.locus.counts,
                        dw->pset->logu,
                        proposal_alpha, dw->pset->prior_alpha, 
                        cumul_aoff,
                        null_sd.sample_points);

    // dw->worker[0]->sample(&null_sd, estimated_mean, dw->final_n_points);
    size_t max_line = 1000;
    while (samples[gs].current != samples[gs].end)
    {
        if (dist_ptr)
        {
            ALLOC_GROW_REMAP(dist_buf->buf, dist_ptr,
                             dist_ptr - dist_buf->buf + max_line, dist_buf->alloc);
            
            dist_ptr = 
                next_distance_quantiles_aux(dw, samples, gs, 
                                            null_sd.sample_points, dist_ptr);
        }
        if (indel_ptr)
        {
            ALLOC_GROW_REMAP(indel_buf->buf, indel_ptr,
                             indel_ptr - indel_buf->buf + max_line,
                             indel_buf->alloc);

            indel_ptr = 
                next_indel_distance_quantiles_aux(dw, samples, gs, indel_ptr);
        }
        if (comp_ptr)
        {
            for (s = 0; s != dw->n_samples; ++s)
            {
                if (! samples[s].is_next) continue;
                if (samples[s].dist_printed)
                {
                    ALLOC_GROW_REMAP(comp_buf->buf, comp_ptr,
                                     comp_ptr - comp_buf->buf + max_line,
                                     comp_buf->alloc);
                    comp_ptr = print_quantiles(dw->pset->comp_quantiles,
                                               dw->pset->n_comp_quantiles,
                                               dw->sample_atts[s].label_string,
                                               &samples[s], comp_ptr);
                    // comp_ptr = dw->worker[s]->print_quantiles(&samples[s], comp_ptr);
                }
            }
        }
        for (s = 0; s != dw->n_samples; ++s)
            if (samples[s].is_next) 
                refresh_locus(&dw->sample_atts[s].nuc_stats,
                              dw->pset->min_quality_score, 
                              &samples[s]);

        gs = init_global_and_next(dw, samples);
    }   

    /* update output buffers */
    if (dist_buf) dist_buf->size = dist_ptr - dist_buf->buf;
    if (comp_buf) comp_buf->size = comp_ptr - comp_buf->buf;
    if (indel_buf) indel_buf->size = indel_ptr - indel_buf->buf;
    
    free(null_sd.sample_points);
    for (s = 0; s != dw->n_samples; ++s)
        free(samples[s].sample_points);

    if (dw->thread_num == 0)
    {
        fprintf(stderr, "Finished processing %s %i\n", 
                samples[gs].locus.reference,
                samples[gs].locus.position);
        fflush(stderr);
    }

    // free(samples);
    delete[] samples;
}


void dist_offload(void *par, const struct managed_buf *bufs)
{
    struct dist_worker_offload_par *ol = (struct dist_worker_offload_par *)par;
    unsigned i = 0;
    if (ol->dist_fh) 
        fwrite(bufs[i].buf, 1, bufs[i].size, ol->dist_fh), i++;

    if (ol->comp_fh)
        fwrite(bufs[i].buf, 1, bufs[i].size, ol->comp_fh), i++;

    if (ol->indel_fh)
        fwrite(bufs[i].buf, 1, bufs[i].size, ol->indel_fh), i++;
}



/*
void print_indel_comparisons(dist_worker_input *wi,
                             locus_sampling *sd,
                             size_t gs)
{

    std::map<unsigned, unsigned>::iterator dit1, e1, dit2, e2;
    CHAR_MAP::iterator iit;

    char contig[100];
    size_t position;
    sscanf(*sd[gs].current, "%s\t%zu", contig, &position);

    for (size_t p = 0; p != wi->n_sample_pairs; ++p)
    {
        size_t s1 = wi->pair_sample1[p], s2 = wi->pair_sample2[p];
        locus_sampling *sd1 = &sd[s1], *sd2 = &sd[s2];

        unsigned max1 = (! sd1->is_next) || sd1->locus.deletions.empty() 
            ? 0 : (*sd1->locus.deletions.rbegin()).first;

        unsigned max2 = (! sd2->is_next) || sd2->locus.deletions.empty()
            ? 0 : (*sd2->locus.deletions.rbegin()).first;

        unsigned max = max1 > max2 ? max1 : max2;

        if (max != 0)
        {
        
            printf("%s\t%s\t%s\t%Zu\t%Zu\t%Zu",
                   wi->worker[s1]->label_string,
                   wi->worker[s2]->label_string,
                   contig,
                   position,
                   sd1->locus.read_depth,
                   sd2->locus.read_depth);



            dit1 = sd1->locus.deletions.begin(), e1 = sd1->locus.deletions.end();
            dit2 = sd2->locus.deletions.begin(), e2 = sd2->locus.deletions.end();

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
