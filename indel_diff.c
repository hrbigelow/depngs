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
indel_event *count_indel_types(struct locus_sampling *sd1,
                               struct locus_sampling *sd2, 
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

        /* now count the non-indel 'events', by definition the
        remaining reads that do not have any indels.  only count this
        as an event if it occurs in at least one of the two pairs.
        read_depth_match are the number of reads with an 'M' base at
        this position.  these may be low-quality basecalls, but the
        call confidence is not relevant here */
        if ((e->count1 = sd1->locus.read_depth_match - nindel1) +
            (e->count2 = sd2->locus.read_depth_match - nindel2) > 0)
        {
            e->seq = NULL;
            ++e;
        }

        *n_counts = e - events;

        return events;
    }    
}


void print_indel_distance_quantiles(const char *contig,
                                    size_t position,
                                    size_t pair_index,
                                    double *dist_quantile_values,
                                    indel_event *events,
                                    size_t n_events,
                                    struct dist_worker_input *dw,
                                    struct managed_buf *mb)
{
    int s1 = sample_pairs.p[pair_index].s1,
        s2 = sample_pairs.p[pair_index].s2;

    unsigned space = (2 * MAX_LABEL_LEN) + 3 + 100 + (10 * MAX_NUM_QUANTILES);
    ALLOC_GROW_TYPED(mb->buf, mb->size + space, mb->alloc);

    mb->size += sprintf(mb->buf + mb->size, 
                        "%s\t%s\t%s\t%Zu", 
                        samples.atts[s1].label,
                        s2 == PSEUDO_SAMPLE ? "REF" : samples.atts[s2].label,
                        contig, position);
    
    for (size_t q = 0; q != n_quantiles; ++q)
        mb->size += sprintf(mb->buf + mb->size, "\t%.4f", dist_quantile_values[q]);

    unsigned indel_space = n_events * (10 + 10);
    ALLOC_GROW_TYPED(mb->buf, mb->size + indel_space, mb->alloc);

    // now print the indel event summary
    indel_event *eb = events, *ee = eb + n_events;
    while (eb != ee) 
    {
        mb->size += sprintf(mb->buf + mb->size, "%c%i",
                            (eb == events ? '\t' : ','), eb->count1);
        eb++;
    }
    eb = events;
    while (eb != ee) 
    {
        mb->size += sprintf(mb->buf + mb->size, "%c%i", 
                            (eb == events ? '\t' : ','), eb->count2);
        eb++;
    }
    eb = events;
    while (eb != ee) 
    {
        unsigned grow = 2 + (eb->seq ? strlen(eb->seq) : 0);
        ALLOC_GROW_TYPED(mb->buf, mb->size + grow, mb->alloc);

        mb->size += sprintf(mb->buf + mb->size, "%c%c%s", 
                            (eb == events ? '\t' : ','), 
                            (eb->seq ? (eb->is_insertion ? '+' : '-') : '@'),
                            (eb->seq ? eb->seq : ""));
        eb++;
    }

    locus_sampling 
        *ls1 = &dw->lslist[s1],
        *ls2 = s2 == PSEUDO_SAMPLE ? &dw->pseudo_sample : &dw->lslist[s2];

    if (worker_options.do_print_pileup)
    {
        unsigned extra_space = 
            ls1->locus.bases_raw.size
            + ls1->locus.quality_codes.size
            + ls2->locus.bases_raw.size
            + ls2->locus.quality_codes.size
            + 50;

        ALLOC_GROW_TYPED(mb->buf, mb->size + extra_space, mb->alloc);
        mb->size += sprintf(mb->buf + mb->size,
                            "\t%Zu\t%s\t%s\t%Zu\t%s\t%s",
                            ls1->locus.read_depth,
                            ls1->locus.bases_raw.buf,
                            ls1->locus.quality_codes.buf,
                            ls2->locus.read_depth,
                            ls2->locus.bases_raw.buf,
                            ls2->locus.quality_codes.buf);
    }
    mb->size += sprintf(mb->buf + mb->size, "\n");
}



// print out all next distance quantiles for indels
void next_indel_distance_quantiles_aux(struct dist_worker_input *dw, 
                                       size_t gs,
                                       struct managed_buf *buf)
{
    /*
      1.  count the indel types
      2.  allocate two buffers
      3.  populate each buffer with dirichlets
      4.  compute the dist quantiles
      5.  deallocate the buffers
      5.  print out suitably filtered distances
    */    
    size_t n_events;

    for (size_t pi = 0; pi != sample_pairs.n; ++pi)
    {
        int s1 = sample_pairs.p[pi].s1, 
            s2 = sample_pairs.p[pi].s2;
        
        struct locus_sampling 
            *ls1 = &dw->lslist[s1], 
            *ls2 = s2 == PSEUDO_SAMPLE ? &dw->pseudo_sample : &dw->lslist[s2];

        if (! (ls1->is_next && ls2->is_next)
            && min_dirichlet_dist > 0) continue;
        
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
        
            size_t bufsize = n_events * max_sample_points;
            double *points1 = new double[bufsize], *p1 = points1, *pe1 = points1 + bufsize;
            double *points2 = new double[bufsize], *p2 = points2;

            while (p1 != pe1) 
            {
                gsl_ran_dirichlet(dw->randgen, n_events, alpha1, p1);
                gsl_ran_dirichlet(dw->randgen, n_events, alpha2, p2);
                p1 += n_events;
                p2 += n_events;
            }

            compute_square_dist(points1, points2, max_sample_points, n_events,
                                dw->square_dist_buf);

            double test_quantile = 1.0 - posterior_confidence, test_quantile_value;
            
            // compute distance quantiles
            compute_marginal_quantiles(dw->square_dist_buf,
                                       max_sample_points,
                                       1, /* one dimensional */
                                       0, /* use the first dimension */
                                       &test_quantile,
                                       1, /* evaluate only one quantile */
                                       &test_quantile_value);

            test_quantile_value = sqrt(test_quantile_value);
            if (test_quantile_value >= min_dirichlet_dist)
            {
                compute_marginal_quantiles(dw->square_dist_buf,
                                           max_sample_points,
                                           1, /* one dimensional */
                                           0, /* use the first dimension */
                                           quantiles,
                                           n_quantiles,
                                           dw->dist_quantile_values);
                unsigned q;
                for (q = 0; q != n_quantiles; ++q)
                    dw->dist_quantile_values[q] = 
                        sqrt(dw->dist_quantile_values[q]) * ONE_OVER_SQRT2;

                print_indel_distance_quantiles(ls1->locus.reference, 
                                               ls1->locus.position, pi, 
                                               dw->dist_quantile_values, 
                                               all_events, n_events, dw, buf);
            }

            delete[] points1;
            delete[] points2;
            delete[] alpha1;
            delete[] alpha2;
        }
        delete[] all_events;
    }
}



#if 0
void print_indel_comparisons(dist_worker_input *dw,
                             locus_sampling *sd,
                             size_t gs)
{

    std::map<unsigned, unsigned>::iterator dit1, e1, dit2, e2;
    CHAR_MAP::iterator iit;

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
                   sd1->locus.reference,
                   sd1->locus.position,
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
