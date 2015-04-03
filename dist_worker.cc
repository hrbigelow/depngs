#include "dist_worker.h"
#include "sampling.h"

#include <gsl/gsl_math.h>
#include <gsl/gsl_randist.h>

extern "C" {
#include "locus.h"
#include "dirichlet_points_gen.h"
}

#define MIN(a, b) ((a) < (b) ? (a) : (b))


void init_sample_attributes(const char *jpd_file,
                            const char *label_string,
                            const char *pileup_file,
                            struct sample_attributes *s)
{
    nucleotide_stats_initialize(jpd_file, &s->nuc_stats);
    if (strlen(label_string) + 1 > sizeof(s->label_string)) 
    {
        fprintf(stderr, "%s: error: sample label string must be "
                "less than %Zu characters\n", __func__,
                sizeof(s->label_string));
        exit(1);
    }
    strcpy(s->label_string, label_string);
    s->fh = fopen(pileup_file, "r");
    if (! s->fh)
    {
        fprintf(stderr, "%s: error: couldn't open pileup input file %s\n",
                __func__, pileup_file);
        exit(1);
    }
}

typedef double COMP_QV[NUM_NUCS][MAX_NUM_QUANTILES];

char *print_basecomp_quantiles(COMP_QV quantile_values,
                               const double *means,
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

    std::multimap<double, size_t, std::greater<double> > dim_to_mean;
    unsigned d, q;
    for (d = 0; d != NUM_NUCS; ++d)
        dim_to_mean.insert(std::make_pair(means[d], d));

    /* calculate mean rank order */
    size_t mean_rank_order[NUM_NUCS];
    std::multimap<double, size_t, std::greater<double> >::iterator dtm;
    d = 0;
    for (dtm = dim_to_mean.begin(); dtm != dim_to_mean.end(); ++dtm)
        mean_rank_order[(*dtm).second] = d++;

    static const char dimension_labels[] = "ACGT";
    for (d = 0; d != NUM_NUCS; ++d)
    {
        out_buffer += 
            sprintf(out_buffer, "%s\t%c\t%Zu\t%10.8f", line_label, 
                    dimension_labels[d], mean_rank_order[d], means[d]);

        for (q = 0; q != n_quantiles; ++q)
            out_buffer += sprintf(out_buffer, "\t%10.8f", quantile_values[d][q]);

        out_buffer += sprintf(out_buffer, "\n");
    }

    return out_buffer;
}


/* computes weighted square euclidean distances between points, using
   the product of their weights as the weight on the final distance. */
void compute_wsq_dist(const POINT *points1,
                      const double *weights1,
                      const POINT *points2,
                      const double *weights2,
                      size_t n_points,
                      double *square_dist_buf,
                      double *weights_buf)
{
    int d;
    double *sd = square_dist_buf, *w = weights_buf;
    size_t np = n_points;
    while (np-- > 0)
    {
        *sd = 0;
        for (d = 0; d != NUM_NUCS; ++d) *sd += gsl_pow_2((*points1)[d] - (*points2)[d]);
        *w++ = *weights1++ * *weights2++;
        sd++;
        points1++;
        points2++;
    }
}


/* unlike previous, this does not require the points to have NUM_NUCS
   components. */
void compute_sq_dist(const double *points1,
                     const double *points2,
                     size_t n_points,
                     size_t n_dims,
                     double *square_dist_buf)
{
    unsigned d;
    double *sd = square_dist_buf;
    size_t np = n_points;
    while (np-- > 0)
    {
        *sd = 0;
        for (d = 0; d != n_dims; ++d) *sd += gsl_pow_2(*points1++ - *points2++);
        ++sd;
    }
}


/* Compute the requested set of distance quantile values from two sets
   of weighted points. */
void compute_dist_wquantiles(double *square_dist_buf,
                             double *weights_buf,
                             size_t n_points,
                             const double *dist_quantiles,
                             size_t n_dist_quantiles,
                             double *dist_quantile_values)
{
    double mean; /* unused */
    compute_marginal_wquantiles(square_dist_buf, weights_buf, n_points, 1, 0,
                                dist_quantiles, n_dist_quantiles,
                                dist_quantile_values,
                                &mean);

    unsigned q;
    for (q = 0; q != n_dist_quantiles; ++q) 
        dist_quantile_values[q] = sqrt(dist_quantile_values[q]);
}


// print out distance quantiles, based on quasi-random pairings of two
// samplings for efficiency, the 'random' pairing is done simply by
// cycling through both sets of sample points, but starting in the
// middle for the second set.  return the next write position.
char *print_distance_quantiles(const char *contig,
                               size_t position,
                               struct dist_worker_input *dw,
                               size_t pair_index,
                               double *dist_quantile_values,
                               locus_sampling *ls,
                               char *out_dist_buf)
{
    size_t s1 = dw->pair_sample1[pair_index], s2 = dw->pair_sample2[pair_index];

    out_dist_buf += sprintf(out_dist_buf, "%s\t%s\t%s\t%Zu", 
                            dw->sample_atts[s1].label_string,
                            dw->sample_atts[s2].label_string,
                            contig, position);
    
    for (size_t q = 0; q != dw->bep.pset->n_dist_quantiles; ++q)
       out_dist_buf += sprintf(out_dist_buf, "\t%7.4f", dist_quantile_values[q]);

    if (dw->print_pileup_fields)
        out_dist_buf += sprintf(out_dist_buf, "\t%Zu\t%s\t%s\t%Zu\t%s\t%s",
                                ls[s1].locus.read_depth,
                                ls[s1].locus.bases_raw.buf,
                                ls[s1].locus.quality_codes.buf,
                                ls[s2].locus.read_depth,
                                ls[s2].locus.bases_raw.buf,
                                ls[s2].locus.quality_codes.buf);

    ls[s1].dist_printed = 1;
    ls[s2].dist_printed = 1;

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



void alloc_pileup_locus(struct locus_sampling *ls,
                        unsigned max_sample_points)
{
    unsigned msp = max_sample_points;
    ls->locus = PileupSummary();
    ls->is_next = false;
    ls->dist_printed = 0;
    ls->distp.pgen = { 
        malloc(sizeof(struct dir_points_par)),
        gen_dirichlet_points_wrapper, 
        malloc(sizeof(struct calc_post_to_dir_par)),
        calc_post_to_dir_ratio
    };
    ((struct dir_points_par *)ls->distp.pgen.point_par)->randgen = gsl_rng_alloc(gsl_rng_taus);
    ls->distp.points = { (POINT *)malloc(sizeof(POINT) * msp), 0, msp };
    ls->distp.weights = { (double *)malloc(sizeof(double) * msp), 0, msp };
    ls->locus_ord = min_pair_ord;
}


void free_pileup_locus(struct locus_sampling *ls)
{
    free((struct dir_points_par *)ls->distp.pgen.point_par);
    free((struct calc_post_to_dir_par *)ls->distp.pgen.weight_par);
    free(ls->distp.points.buf);
    free(ls->distp.weights.buf);
    gsl_rng_free(((struct dir_points_par *)ls->distp.pgen.point_par)->randgen);
}


/* update 'ls' fields to be consistent with sd->current */
void update_pileup_locus(const struct nucleotide_stats *stats,
                         size_t min_quality_score,
                         struct dist_worker_input *dw,
                         locus_sampling *ls)
{
    ls->locus.load_line(ls->current);
    ls->locus_ord = init_locus(ls->current);
    ls->locus.parse(min_quality_score);
    ls->distp.points.size = 0;
    ls->distp.weights.size = 0;
    nucleotide_stats_pack(stats, &ls->locus.counts);
    struct calc_post_to_dir_par *wp = 
        (struct calc_post_to_dir_par *)ls->distp.pgen.weight_par;

    wp->post_counts = &ls->locus.counts;
    (void)alphas_from_counts(wp->post_counts,
                             dw->bep.pset->prior_alpha,
                             wp->proposal_alpha);
    struct dir_points_par *pp = 
        (struct dir_points_par *)ls->distp.pgen.point_par;

    memcpy(pp->alpha, wp->proposal_alpha, sizeof(pp->alpha));
}


/* advance ls->current, then initialize if there is another locus */
void refresh_locus(const struct nucleotide_stats *stats,
                   size_t min_quality_score,
                   struct dist_worker_input *dw,
                   locus_sampling *ls)
{
    assert(ls->current != ls->end);
    if ((ls->current = strchr(ls->current, '\n') + 1) == ls->end) 
        ls->is_next = false;

    else
    {
        update_pileup_locus(stats, min_quality_score, dw, ls);
        ls->dist_printed = 0;
    }
}


char *print_indel_distance_quantiles(const char *contig,
                                     size_t position,
                                     dist_worker_input *dw,
                                     size_t pair_index,
                                     double *dist_quantile_values,
                                     indel_event *events,
                                     size_t n_events,
                                     locus_sampling *sd,
                                     char *out_dist_buf)
{
    size_t s1 = dw->pair_sample1[pair_index], s2 = dw->pair_sample2[pair_index];

    out_dist_buf += sprintf(out_dist_buf, "%s\t%s\t%s\t%Zu", 
                            dw->sample_atts[s1].label_string,
                            dw->sample_atts[s2].label_string,
                            contig, position);
    
    for (size_t q = 0; q != dw->bep.pset->n_dist_quantiles; ++q)
        out_dist_buf += sprintf(out_dist_buf, "\t%.4f", dist_quantile_values[q]);

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
char *next_distance_quantiles_aux(dist_worker_input *dw, 
                                  locus_sampling *lslist,
                                  size_t gs,
                                  char *out_buf)
{
    if (! out_buf) return out_buf;

    locus_sampling *lsp[2];
    size_t pi, i;
    unsigned counts1[NUM_NUCS], counts2[NUM_NUCS];
    enum fuzzy_state diff_state = AMBIGUOUS;


    for (pi = 0; pi != dw->n_sample_pairs; ++pi)
    {
        dw->bep.dist[0] = &lslist[dw->pair_sample1[pi]].distp;
        dw->bep.dist[1] = &lslist[dw->pair_sample2[pi]].distp;
        lsp[0] = &lslist[dw->pair_sample1[pi]];
        lsp[1] = &lslist[dw->pair_sample2[pi]];

        /* TEMPORARY */
        print_beb_bounds(&dw->bep);
        exit(1);

        for (i = 0; i != NUM_NUCS; ++i)
        {
            counts1[i] = lsp[0]->locus.base_counts[i];
            counts2[i] = lsp[1]->locus.base_counts[i];
        }

        if (! (lsp[0]->is_next && lsp[1]->is_next))
            diff_state = AMBIGUOUS;

        else
            diff_state = cached_dirichlet_diff(counts1, counts2, &dw->bep);
        
        if (diff_state == CHANGED)
        {
            /* Finish sampling and do full distance marginal estimation */
            for (i = 0; i != 2; ++i)
            {
                struct distrib_points *dst = dw->bep.dist[i];
                /* Generate missing points */
                POINT *p,
                    *pb = dst->points.buf + dst->points.size,
                    *pe = dst->points.buf + dw->bep.pset->max_sample_points;
                for (p = pb; p != pe; p += GEN_POINTS_BATCH)
                    dst->pgen.gen_point(dst->pgen.point_par, p);

                dst->points.size = dw->bep.pset->max_sample_points;
                
                /* Generate missing weights */
                double *w,
                    *wb = dst->weights.buf + dst->weights.size,
                    *we = dst->weights.buf + dw->bep.pset->max_sample_points;
                for (w = wb, p = dst->points.buf + dst->weights.size; w != we; 
                     w += GEN_POINTS_BATCH, p += GEN_POINTS_BATCH)
                    dst->pgen.weight(p, dst->pgen.weight_par, w);
            }
            /* Compute weighted square distances */
            compute_wsq_dist(dw->bep.dist[0]->points.buf, 
                             dw->bep.dist[0]->weights.buf,
                             dw->bep.dist[1]->points.buf, 
                             dw->bep.dist[1]->weights.buf,
                             dw->bep.pset->max_sample_points,
                             dw->square_dist_buf, dw->weights_buf);

            double test_quantile = 1.0 - dw->bep.pset->post_confidence, test_quantile_value;

            /* Compute the test distance quantile */
            compute_dist_wquantiles(dw->square_dist_buf,
                                    dw->weights_buf,
                                    dw->bep.pset->max_sample_points,
                                    &test_quantile,
                                    1,
                                    &test_quantile_value);

            test_quantile_value = sqrt(test_quantile_value);
            if (test_quantile_value > dw->bep.pset->min_dist)
            {
                compute_dist_wquantiles(dw->square_dist_buf,
                                        dw->weights_buf,
                                        dw->bep.pset->max_sample_points,
                                        dw->bep.pset->dist_quantiles,
                                        dw->bep.pset->n_dist_quantiles,
                                        dw->dist_quantile_values);            

                unsigned q;
                for (q = 0; q != dw->bep.pset->n_dist_quantiles; ++q)
                    dw->dist_quantile_values[q] = sqrt(dw->dist_quantile_values[q]);

                out_buf = 
                    print_distance_quantiles(lslist[gs].locus.reference, 
                                             lslist[gs].locus.position, dw, 
                                             pi, dw->dist_quantile_values, 
                                             lslist, out_buf);
            }
        }
        
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
char *next_indel_distance_quantiles_aux(dist_worker_input *dw, 
                                        locus_sampling *lslist,
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
    
    sscanf(lslist[gs].current, "%s\t%zu", contig, &position);

    for (size_t pi = 0; pi != dw->n_sample_pairs; ++pi)
    {
        size_t s1 = dw->pair_sample1[pi], s2 = dw->pair_sample2[pi];
        locus_sampling *ls1 = &lslist[s1], *ls2 = &lslist[s2];

        if (! (ls1->is_next && ls2->is_next) && dw->bep.pset->min_dist > 0) continue;
        
        indel_event *all_events = count_indel_types(ls1, ls2, &n_events);

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
        
            size_t bufsize = n_events * dw->bep.pset->max_sample_points;
            double *points1 = new double[bufsize], *p1 = points1, *pe1 = points1 + bufsize;
            double *points2 = new double[bufsize], *p2 = points2;

            while (p1 != pe1) 
            {
                gsl_ran_dirichlet(dw->randgen, n_events, alpha1, p1);
                gsl_ran_dirichlet(dw->randgen, n_events, alpha2, p2);
                p1 += n_events;
                p2 += n_events;
            }

            compute_sq_dist(points1, points2, dw->bep.pset->max_sample_points, n_events,
                            dw->square_dist_buf);

            double test_quantile = 1.0 - dw->bep.pset->post_confidence, test_quantile_value;
            
            // compute distance quantiles
            compute_marginal_quantiles(dw->square_dist_buf,
                                       dw->bep.pset->max_sample_points,
                                       1, /* one dimensional */
                                       0, /* use the first dimension */
                                       &test_quantile,
                                       1, /* evaluate only one quantile */
                                       &test_quantile_value);

            test_quantile_value = sqrt(test_quantile_value);
            if (test_quantile_value >= dw->bep.pset->min_dist)
            {
                compute_marginal_quantiles(dw->square_dist_buf,
                                           dw->bep.pset->max_sample_points,
                                           1, /* one dimensional */
                                           0, /* use the first dimension */
                                           dw->bep.pset->dist_quantiles,
                                           dw->bep.pset->n_dist_quantiles,
                                           dw->dist_quantile_values);
                unsigned q;
                for (q = 0; q != dw->bep.pset->n_dist_quantiles; ++q)
                    dw->dist_quantile_values[q] = sqrt(dw->dist_quantile_values[q]);

                out_buf = 
                    print_indel_distance_quantiles(contig, position, dw, pi, 
                                                   dw->dist_quantile_values, 
                                                   all_events, n_events, lslist, out_buf);
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
size_t init_global_and_next(dist_worker_input *dw, locus_sampling *samples)
{
    size_t gs = 0;
    for (size_t s = 0; s != dw->n_samples; ++s)
    {
        if (samples[s].current != samples[s].end
            && (samples[gs].current == samples[gs].end
                || cmp_pair_ordering(&samples[s].locus_ord, 
                                     &samples[gs].locus_ord) < 0)
            )
            gs = s;
    }
    bool global_at_end = (samples[gs].current == samples[gs].end);
    for (size_t s = 0; s != dw->n_samples; ++s)
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
    struct locus_sampling *lslist = new struct locus_sampling[dw->n_samples];

    size_t gs = 0, s;

    for (s = 0; s != dw->n_samples; ++s)
    {    
        alloc_pileup_locus(&lslist[s], dw->bep.pset->max_sample_points);
        lslist[s].current = in_bufs[s].buf;
        lslist[s].end = in_bufs[s].buf + in_bufs[s].size;
    }

    
    // before main loop, initialize loci and samples
    for (s = 0; s != dw->n_samples; ++s)
    {
        if (lslist[s].current == lslist[s].end) continue;
        update_pileup_locus(&dw->sample_atts[s].nuc_stats, 
                            dw->bep.pset->min_quality_score, 
                            dw,
                            &lslist[s]);
    }

    gs = init_global_and_next(dw, lslist);

    // main loop for computing pairwise distances
    // though slightly inefficient, just construct null_points here
    // !!! here, use worker[0] as a proxy.  this is a design flaw
    // owing to the fact that there are multiple workers used, all with the same
    // basic settings
    locus_sampling null_sd;
    alloc_pileup_locus(&null_sd, dw->bep.pset->max_sample_points);

    char null_pileup[] = "chr1\t1\tA\t1\t\n";
    null_sd.current = null_pileup;
    null_sd.end = null_pileup + sizeof(null_pileup);
        
    update_pileup_locus(&dw->sample_atts[0].nuc_stats, 
                        dw->bep.pset->min_quality_score,
                        dw, &null_sd);

    COMP_QV comp_quantile_values;
    double comp_means[NUM_NUCS];
    
    size_t max_line = 1000;
    while (lslist[gs].current != lslist[gs].end)
    {
        if (dist_ptr)
        {
            ALLOC_GROW_REMAP(dist_buf->buf, dist_ptr,
                             dist_ptr - dist_buf->buf + max_line, dist_buf->alloc);
            
            dist_ptr = next_distance_quantiles_aux(dw, lslist, gs, dist_ptr);
        }
        if (indel_ptr)
        {
            ALLOC_GROW_REMAP(indel_buf->buf, indel_ptr,
                             indel_ptr - indel_buf->buf + max_line,
                             indel_buf->alloc);

            indel_ptr = next_indel_distance_quantiles_aux(dw, lslist, gs, indel_ptr);
        }
        if (comp_ptr)
        {
            for (s = 0; s != dw->n_samples; ++s)
            {
                if (! lslist[s].is_next) continue;
                if (lslist[s].dist_printed)
                {
                    ALLOC_GROW_REMAP(comp_buf->buf, comp_ptr,
                                     comp_ptr - comp_buf->buf + max_line,
                                     comp_buf->alloc);

                    unsigned d;
                    for (d = 0; d != NUM_NUCS; ++d)
                        compute_marginal_wquantiles((double *)lslist[s].distp.points.buf,
                                                    lslist[s].distp.weights.buf,
                                                    dw->bep.pset->max_sample_points,
                                                    NUM_NUCS,
                                                    d,
                                                    dw->bep.pset->comp_quantiles,
                                                    dw->bep.pset->n_comp_quantiles,
                                                    comp_quantile_values[d],
                                                    &comp_means[d]);
                    
                    comp_ptr = print_basecomp_quantiles(comp_quantile_values,
                                                        comp_means,
                                                        dw->bep.pset->n_comp_quantiles,
                                                        dw->sample_atts[s].label_string,
                                                        &lslist[s],
                                                        comp_ptr);
                }
            }
        }
        for (s = 0; s != dw->n_samples; ++s)
            if (lslist[s].is_next) 
                refresh_locus(&dw->sample_atts[s].nuc_stats,
                              dw->bep.pset->min_quality_score, 
                              dw,
                              &lslist[s]);

        gs = init_global_and_next(dw, lslist);
    }   

    /* update output buffers */
    if (dist_buf) dist_buf->size = dist_ptr - dist_buf->buf;
    if (comp_buf) comp_buf->size = comp_ptr - comp_buf->buf;
    if (indel_buf) indel_buf->size = indel_ptr - indel_buf->buf;
    
    if (dw->thread_num == 0)
    {
        fprintf(stderr, "Finished processing %s %i\n", 
                lslist[gs].locus.reference,
                lslist[gs].locus.position);
        fflush(stderr);
    }

    for (s = 0; s != dw->n_samples; ++s) free_pileup_locus(&lslist[s]);

    free_pileup_locus(&null_sd);

    delete[] lslist;
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





/* ***************************************************************** */



#if 0
void print_indel_comparisons(dist_worker_input *dw,
                             locus_sampling *sd,
                             size_t gs)
{

    std::map<unsigned, unsigned>::iterator dit1, e1, dit2, e2;
    CHAR_MAP::iterator iit;

    char contig[100];
    size_t position;
    sscanf(*sd[gs].current, "%s\t%zu", contig, &position);

    for (size_t p = 0; p != dw->n_sample_pairs; ++p)
    {
        size_t s1 = dw->pair_sample1[p], s2 = dw->pair_sample2[p];
        locus_sampling *sd1 = &sd[s1], *sd2 = &sd[s2];

        unsigned max1 = (! sd1->is_next) || sd1->locus.deletions.empty() 
            ? 0 : (*sd1->locus.deletions.rbegin()).first;

        unsigned max2 = (! sd2->is_next) || sd2->locus.deletions.empty()
            ? 0 : (*sd2->locus.deletions.rbegin()).first;

        unsigned max = max1 > max2 ? max1 : max2;

        if (max != 0)
        {
        
            printf("%s\t%s\t%s\t%Zu\t%Zu\t%Zu",
                   dw->worker[s1]->label_string,
                   dw->worker[s2]->label_string,
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
#endif
