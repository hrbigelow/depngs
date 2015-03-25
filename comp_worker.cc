#include "comp_worker.h"

#include <gsl/gsl_sf_exp.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

#include "sampling.h"
#include "pileup_tools.h"
#include "defs.h"

#include "yepMath.h"
#include "yepCore.h"

extern "C" {
#include "nucleotide_stats.h"
#include "cache.h"
#include "locus.h"
#include "binomial_est.h"
#include "metropolis_sampling.h"
}


#if 0
/* returns next write position after writing to out_buffer */
char *process_line_comp_binom(struct comp_worker_input *cw,
                              struct locus_sampling *ls,
                              struct points_buf *ref_points,
                              struct weights_buf *ref_weights,
                              char *out_buffer)
{
    /* use this as the dummy reference point generator */
    struct points_gen ref_gen = {
        &ls->locus.reference_base,
        gen_reference_points_wrapper,
        NULL,
        calc_dummy_ratio
    };

    ref_points->size = ref_weights->size = 0;

    enum fuzzy_state sampling_state =
        binomial_quantile_est(cw->pset.max_sample_points,
                              cw->pset.min_dist,
                              cw->pset.post_confidence, 
                              cw->pset.beta_confidence,
                              ls->pgen,
                              &ls->points,
                              &ls->weights,
                              ref_gen,
                              ref_points,
                              ref_weights,
                              GEN_POINTS_BATCH);
                              
    if (sampling_state == CHANGED)
    {
        /* Finish the sampling and do a full marginal estimation */
        POINT *p,
            *pb = ls->points.buf + ls->points.size,
            *pe = ls->points.buf + cw->pset.max_sample_points;
        for (p = pb; p != pe; p += GEN_POINTS_BATCH)
            ls->pgen.gen_point(ls->pgen.gen_point_par, p);
        ls->points.size = cw->pset.max_sample_points;

        /* Generate all missing weights.  This is the expensive part */
        double *w,
            *wb = ls->weights.buf + ls->weights.size,
            *we = ls->weights.buf + cw->pset.max_sample_points;
        for (w = wb, p = ls->points.buf + ls->weights.size; w != we; ++w, ++p)
            ls->pgen.weight(p, ls->pgen.weight_par, w);

        int ref_dim = base_to_index(ls->locus.reference_base);
        double test_quantile = 1.0 - cw->pset.post_confidence, test_quantile_value;

        compute_marginal_wquantiles(ls->points.buf,
                                    ls->weights.buf,
                                    ls->points.size,
                                    ref_dim,
                                    &test_quantile,
                                    1,
                                    &test_quantile_value);

        if (test_quantile_value > cw->pset.min_dist)
        {
            double quantile_values[NUM_NUCS][MAX_NUM_QUANTILES];
            int d;
            for (d = 0; d != NUM_NUCS; ++d)
                compute_marginal_wquantiles(ls->points.buf,
                                            ls->weights.buf,
                                            ls->points.size,
                                            d,
                                            cw->quantiles,
                                            cw->n_quantiles,
                                            quantile_values[d]);
            
            out_buffer = print_quantiles(cw->quantiles,
                                         cw->n_quantiles,
                                         cw->sample_atts.label_string,
                                         ls, out_buffer);
        }
    }
    return out_buffer;
}
#endif


#if 0
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

    struct points_buf ref_points = { NULL, 0, cw->pset.max_sample_points };
    struct points_buf loc_points = { NULL, 0, cw->pset.max_sample_points };
    
    struct weights_buf ref_weights = { NULL, 0, cw->pset.max_sample_points };
    struct weights_buf loc_weights = { NULL, 0, cw->pset.max_sample_points };
        
    ALLOC_GROW_TYPED(ref_points.buf, ref_points.size, ref_points.alloc);
    ALLOC_GROW_TYPED(ref_weights.buf, ref_weights.size, ref_weights.alloc);
    ALLOC_GROW_TYPED(loc_points.buf, loc_points.size, loc_points.alloc);
    ALLOC_GROW_TYPED(loc_weights.buf, loc_weights.size, loc_weights.alloc);
    
    struct calc_post_to_dir_par post_dir_par;
    struct dir_points_par dir_par;
    dir_par.randgen = gsl_rng_alloc(gsl_rng_taus);
    
    struct locus_sampling ls = {
        PileupSummary(),
        false,
        0,
        { &dir_par, gen_dirichlet_points_wrapper, &post_dir_par, calc_post_to_dir_ratio },
        loc_points,
        loc_weights,
        in_buf->buf,
        in_buf->buf + in_buf->size,
        min_pair_ord
    };
    post_dir_par.post_counts = &ls.locus.counts;
    
    init_pileup_locus(&cw->sample_atts.nuc_stats, cw->pset.min_quality_score, &ls);
    
    while (ls.current != ls.end)
    {
        ALLOC_GROW_REMAP(out_buf->buf, write_ptr, 
                         write_ptr - out_buf->buf + max_line, out_buf->alloc);
        
        write_ptr = process_line_comp_binom(cw, &ls, &ref_points, &ref_weights, write_ptr);
        
        refresh_locus(&cw->sample_atts.nuc_stats,
                      cw->pset.min_quality_score, &ls);
    }
    out_buf->size = write_ptr - out_buf->buf;
    
    gsl_rng_free(dir_par.randgen);
    free(ls.points.buf);
    free(ls.weights.buf);
    free(ref_points.buf);
    free(ref_weights.buf);
}
#endif



/* returns next write position after writing to out_buffer
   sample_points_buf is a scratch space for writing sample points must
   be this->final_n_points * 4 size */
#if 0
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
#endif


