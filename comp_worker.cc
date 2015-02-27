#include "comp_worker.h"

#include <gsl/gsl_sf_exp.h>

#include "sampling.h"
#include "pileup_tools.h"
#include "defs.h"

extern "C" {
#include "nucleotide_stats.h"
#include "cache.h"
#include "locus.h"
#include "nucleotide_stats.h"
}

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



/* initialize the locus in 'ls' defined by sd->current */
void init_pileup_locus(const struct nucleotide_stats *stats,
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
    /*
    const struct cpd_count
        *pc = ls->locus.counts.stats,
        *pce = pc + ls->locus.counts.num_data;
    while (pc != pce)
    {
        fprintf(stderr, "{ %g, %g, %g, %g, %lu },\n",
        pc->cpd[0], pc->cpd[1], pc->cpd[2], pc->cpd[3], pc->ct);
        ++pc;
    }
    */

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
    struct eval_counts eval;
    size_t cumul_aoff = 
        tune_proposal(&ls->locus.counts,
                      &cw->pset, proposal_alpha, estimated_mean,
                      ls->sample_points, &eval);
    // this->tune(&ls, initial_point);
    metropolis_sampling(0, cw->pset.final_n_points, &ls->locus.counts,
                        cw->pset.logu, proposal_alpha, 
                        cw->pset.prior_alpha, cumul_aoff,
                        ls->sample_points, &eval);

    ls->n_sample_points = cw->pset.final_n_points;

    // this->sample(&ls, initial_point, this->s.final_n_points);

    // we always print a base composition estimate for loci with
    // ref=N.  So, we only do the filtering test if ref!=N
    if (ref_ind >= 0)
        compute_marginal_quantiles(ls->sample_points,
                                   ls->n_sample_points,
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
    struct comp_worker_input *cw = (struct comp_worker_input *)par;

    char *write_ptr = out_buf->buf;
    size_t max_line = 1000;
    struct locus_sampling ls = {
        PileupSummary(),
        false,
        0,
        (double *)malloc(cw->pset.final_n_points * NUM_NUCS * sizeof(double)),
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
