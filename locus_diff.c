#define __STDC_LIMIT_MACROS
#include <stdint.h>

#include "locus_diff.h"

#include "wquantile.h"
#include "sampling.h"

#include <gsl/gsl_math.h>
#include <gsl/gsl_randist.h>
#include <time.h>
#include <pthread.h>
#include <assert.h>

#include "bam_sample_info.h"
#include "chunk_strategy.h"
#include "geometry.h"
#include "fasta.h"
#include "timer.h"
#include "dir_cache.h"
#include "ksort.h"

#define MIN(a, b) ((a) < (b) ? (a) : (b))

static struct {
    struct bam_scanner_info *reader_buf;
    void **reader_pars;
    unsigned n_threads;
    struct locus_diff_offload_par offload_par;
    const char *fasta_file;
} thread_params;

static struct locus_diff_params g_ld_par;

static __thread struct locus_diff_input tls_dw;


#define NUM_EXTRA_BUFS_PER_THREAD 50

void
locus_diff_thread_init();

void
locus_diff_thread_free();


struct work_unit {
    unsigned n_threads, nth;
};


/* initialize all bam_stats for sample nth and thread 0, and duplicate
   them for all other threads. */
void *
stats_init_worker(void *arg)
{
    struct work_unit *wu = arg;
    unsigned t = wu->nth;

    /* once-per-thread initialization */
    thread_params.reader_buf[t] = (struct bam_scanner_info){
        .m = malloc(bam_samples.n * sizeof(struct bam_stats)),
        .n = bam_samples.n,
        .do_print_progress = 1
    };
    thread_params.reader_pars[t] = &thread_params.reader_buf[t];
    
    unsigned s, dt;
    for (s = t; s < bam_samples.n; s += wu->n_threads) {
        struct bam_stats *bstats =
            &thread_params.reader_buf[0].m[s];
        bam_stats_init(bam_samples.m[s].bam_file, bstats);
        for (dt = 1; dt != wu->n_threads; ++dt)
            thread_params.reader_buf[dt].m[s] =
                bam_stats_dup(bstats, bam_samples.m[s].bam_file);
    }

    return NULL;
}


struct thread_queue *
locus_diff_init(const char *samples_file, const char *sample_pairs_file, 
                const char *locus_range_file, const char *fasta_file,
                unsigned n_threads, unsigned n_max_reading, unsigned long max_input_mem,
                struct locus_diff_params ld_par,
                struct binomial_est_params be_par,
                struct dir_cache_params dc_par,
                struct bam_filter_params bf_par,
                FILE *dist_fh, FILE *comp_fh, FILE *indel_fh)
{
    g_ld_par = ld_par;
    if (g_ld_par.post_confidence < 0.8
        || g_ld_par.post_confidence > 0.999999) {
        fprintf(stderr,
                "Error: locus_diff_init: posterior confidence of %g is a bad value\n",
                g_ld_par.post_confidence);
        exit(1);
    }
    if (g_ld_par.min_dirichlet_dist <= 0 || g_ld_par.min_dirichlet_dist >= 1) {
        fprintf(stderr,
                "Error: locus_diff_init: min_dirichlet_dist of %g is a bad value.\n"
                "Should be in [0, 1]\n", g_ld_par.min_dirichlet_dist);
        exit(1);
    }
    /* */
    unsigned long 
        cs_bytes_zone2 = 1e8, 
        cs_bytes_zone3 = 1e7;

    bam_sample_info_init(samples_file, sample_pairs_file);
    chunk_strategy_init(bam_samples.n, n_threads, 
                        locus_range_file, fasta_file,
                        cs_bytes_zone2, cs_bytes_zone3);

    thread_params.n_threads = n_threads;
    thread_params.reader_buf = malloc(n_threads * sizeof(struct bam_scanner_info));
    thread_params.reader_pars = malloc(n_threads * sizeof(void *));
    
    /* initialize bam_stats */
    fprintf(stdout, "%s: Starting reading BAM indices\n", timer_progress());   
    unsigned t;
    struct work_unit *stats_init_input = malloc(n_threads * sizeof(struct work_unit));
    pthread_t *bam_init_th = malloc(n_threads * sizeof(pthread_t));
    for (t = 0; t != n_threads; ++t) {
        stats_init_input[t] = (struct work_unit){ n_threads, t };
        (void)pthread_create(&bam_init_th[t], NULL,
                             stats_init_worker, &stats_init_input[t]);
    }
    for (t = 0; t != n_threads; ++t)
        pthread_join(bam_init_th[t], NULL);
        
    free(bam_init_th);
    free(stats_init_input);
    fprintf(stdout, "%s: Finished reading BAM indices.\n", timer_progress());
    
    thread_params.fasta_file = fasta_file;

    dir_cache_init(be_par, dc_par, bf_par,
                   thread_params.reader_buf,
                   n_max_reading, max_input_mem, n_threads);

    thread_params.offload_par = 
        (struct locus_diff_offload_par){ dist_fh, comp_fh, indel_fh };

    /* To avoid a stall, n_extra / n_threads should be greater than
       Max(work chunk time) / Avg(work chunk time). */
    size_t n_extra = n_threads * NUM_EXTRA_BUFS_PER_THREAD;
    size_t n_output_files = (dist_fh ? 1 : 0) + (comp_fh ? 1 : 0) + (indel_fh ? 1 : 0);

    thread_queue_reader_t reader = { bam_reader, bam_scanner };

    /* we do not want to skip empty loci, because we need to traverse
       these in order to get statistics for missing data */
    unsigned skip_empty_loci = 0;
    batch_pileup_init(bf_par, skip_empty_loci, dc_par.pseudo_depth);

    /* now turn on progress messages */
    for (t = 0; t != n_threads; ++t)
        thread_params.reader_buf[t].do_print_progress = 1;
    
    struct thread_queue *tqueue =
        thread_queue_init(reader, 
                          thread_params.reader_pars,
                          locus_diff_worker,
                          locus_diff_offload, 
                          &thread_params.offload_par,
                          locus_diff_thread_init,
                          locus_diff_thread_free,
                          n_threads, n_extra, n_max_reading, bam_samples.n,
                          n_output_files, max_input_mem);
    return tqueue;
}


void
locus_diff_free(struct thread_queue *tq)
{
    dir_cache_free();
    batch_pileup_free();
    chunk_strategy_free();
    bam_sample_info_free();

    unsigned t, s;
    for (t = 0; t != thread_params.n_threads; ++t) {
        for (s = 0; s != bam_samples.n; ++s)
            bam_stats_free(&thread_params.reader_buf[t].m[s]);
        free(thread_params.reader_buf[t].m);
    }

    free(thread_params.reader_buf);
    free(thread_params.reader_pars);
    thread_queue_free(tq);
}



/* called just after this thread is created */
void
locus_diff_thread_init()
{
    unsigned msp = g_ld_par.max_sample_points;
    tls_dw.randgen = gsl_rng_alloc(gsl_rng_taus);
    alloc_locus_data(&tls_dw.pseudo_sample);

    /* pseudo-sample weights are filled once, and only here. */
    unsigned i;
    for (i = 0; i != msp; ++i)
        tls_dw.pseudo_sample.dist.weights[i] = 1.0;
    
    tls_dw.ldat = malloc(bam_samples.n * sizeof(struct locus_data));
    unsigned s;
    for (s = 0; s != bam_samples.n; ++s)
        alloc_locus_data(&tls_dw.ldat[s]);

    tls_dw.pair_stats = calloc(bam_sample_pairs.n, sizeof(struct pair_dist_stats));
    tls_dw.square_dist_buf = malloc(sizeof(double) * msp);
    tls_dw.weights_buf = malloc(sizeof(double) * msp);

    batch_pileup_thread_init(bam_samples.n, 
                             thread_params.fasta_file);
    dir_points_thread_init();
}


/* when a worker thread exits, it must inform the rest of the program
   that it will not be modifying shared data anymore. */
void
locus_diff_thread_free()
{
    dir_points_thread_free();
    gsl_rng_free(tls_dw.randgen);

    free_locus_data(&tls_dw.pseudo_sample);
    unsigned s;
    for (s = 0; s != bam_samples.n; ++s)
        free_locus_data(&tls_dw.ldat[s]);
    free(tls_dw.ldat);

    free(tls_dw.pair_stats);
    free(tls_dw.square_dist_buf);
    free(tls_dw.weights_buf);

    batch_pileup_thread_free();
}




struct dim_mean {
    unsigned dim;
    double mean;
};

int
less_mean(const void *ap, const void *bp)
{
    const struct dim_mean *a = ap, *b = bp;
    return a->mean < b->mean;
}
 
void
print_basecomp_quantiles(comp_quantile_vals_t quantile_values,
                         const double *means,
                         size_t n_quantiles,
                         const char *label_string,
                         struct pileup_locus_info *pli,
                         struct pileup_data *pdat,
                         struct managed_buf *mb)
{
    char line_label[2048];

    sprintf(line_label,
            "%s\t%s\t%i\t%c\t%u\t%u\t%u",
            label_string, 
            pli->refname,
            pli->pos + 1,
            pli->refbase,
            pdat->n_match_hi_q,
            pdat->n_match_lo_q,
            pdat->n_indel);

    struct dim_mean dim_to_mean[NUM_NUCS];

    unsigned d;
    for (d = 0; d != NUM_NUCS; ++d)
        dim_to_mean[d] = (struct dim_mean){ d, means[d] };

    qsort(dim_to_mean, NUM_NUCS, sizeof(dim_to_mean[0]), less_mean);

    /* calculate mean rank order */
    size_t mean_rank_order[NUM_NUCS];
    for (d = 0; d != NUM_NUCS; ++d)
        mean_rank_order[dim_to_mean[d].dim] = d;

    unsigned grow = NUM_NUCS * (sizeof(line_label) + 30 + (11 * MAX_NUM_QUANTILES));
    ALLOC_GROW(mb->buf, mb->size + grow, mb->alloc);

    static const char dimension_labels[] = "ACGT";
    for (d = 0; d != NUM_NUCS; ++d) {
        mb->size += sprintf(mb->buf + mb->size,
                            "%s\t%c\t%Zu\t%10.8f", line_label, 
                            dimension_labels[d], mean_rank_order[d], means[d]);
        
        unsigned q;
        for (q = 0; q != n_quantiles; ++q)
            mb->size += sprintf(mb->buf + mb->size, "\t%10.8f", quantile_values[d][q]);

        mb->size += sprintf(mb->buf + mb->size, "\n");
    }
}


/* print out distance quantiles. use a pseudo-sample for the second
   one if its pair index is equal to REFERENCE_SAMPLE */
void
print_distance_quantiles(const char *contig,
                         size_t position,
                         char ref_base,
                         unsigned pair_index,
                         double *dist_quantile_values,
                         struct managed_buf *buf)
{
    unsigned s[] = { bam_sample_pairs.m[pair_index].s1,
                     bam_sample_pairs.m[pair_index].s2 };

    unsigned space = (2 * MAX_LABEL_LEN) + 3 + 100 + (10 * MAX_NUM_QUANTILES);

    ALLOC_GROW(buf->buf, buf->size + space, buf->alloc);

    buf->size +=
        sprintf(buf->buf + buf->size, 
                "%s\t%s\t%s\t%Zu\t%c", 
                bam_samples.m[s[0]].label,
                s[1] == REFERENCE_SAMPLE ? "REF" : bam_samples.m[s[1]].label,
                contig,
                position,
                ref_base);
    
    unsigned q;
    for (q = 0; q != g_ld_par.n_quantiles; ++q)
        buf->size += sprintf(buf->buf + buf->size, "\t%7.4f",
                             dist_quantile_values[q]);

    if (g_ld_par.do_print_pileup) {
        struct locus_data *ld[] = { 
            &tls_dw.ldat[s[0]],
            s[1] == REFERENCE_SAMPLE ? &tls_dw.pseudo_sample : &tls_dw.ldat[s[1]]
        };

        unsigned i;
        for (i = 0; i != 2; ++i)
            if (! ld[i]->init.sample_data) {
                pileup_current_data(s[i], &ld[i]->sample_data);
                ld[i]->init.sample_data = 1;
            }
        
        struct pileup_data *pd[] = { &ld[0]->sample_data, &ld[1]->sample_data };
        unsigned extra_space = 
            pd[0]->calls.size + pd[0]->quals.size
            + pd[1]->calls.size + pd[1]->quals.size
            + 50;
        
        ALLOC_GROW(buf->buf, buf->size + extra_space, buf->alloc);
        
        char *out = buf->buf + buf->size;
        for (i = 0; i != 2; ++i) {
            *out++ = '\t';
            out += sprintf(out, "%u", pd[i]->n_match_hi_q);
            *out++ = '\t';
            strncpy(out, pd[i]->calls.buf, pd[i]->calls.size);
            out += pd[i]->calls.size;
            *out++ = '\t';
            strncpy(out, pd[i]->quals.buf, pd[i]->quals.size);
            out += pd[i]->quals.size;
        }
        *out++ = '\n';
        buf->size = out - buf->buf;
    }
}


struct pos_pair_key {
    struct contig_pos pos;
    struct dir_cache_pair_key pk;
    enum fuzzy_state state;
};


struct pos_pair_iter {
    struct pos_pair_key *buf, *cur;
    unsigned size, alloc;
};


static unsigned
is_pseudo_sample(struct dir_points *dp)
{
    return (dir_cache_is_pseudo_alpha(dp->perm_alpha) != -1);
}


#define ONE_OVER_SQRT2 0.70710678118654752440

/* for each sample pair, calculate whether the loci differ above the
   given level and confidence. if a pair differs, set the
   confirmed_changed flag for each sample in the pair. if out_buf is
   not NULL, print out distance quantiles.  also generates sample
   points for each sample as needed, both for the preliminary test and
   more points for the final test */
void
distance_quantiles_aux(struct managed_buf *out_buf,
                       struct pos_pair_iter *dist_class)
{
    enum fuzzy_state diff_state = AMBIGUOUS;
    unsigned cache_hit;
    struct locus_data *ld[2];
    unsigned pi, i;
    struct contig_pos pos = pileup_current_pos();
    
    for (pi = 0; pi != bam_sample_pairs.n; ++pi) {
        struct pos_pair_iter *it = &dist_class[pi];
        assert(cmp_contig_pos(pos, it->cur->pos) == 0);
        
        unsigned sp[] = { bam_sample_pairs.m[pi].s1, bam_sample_pairs.m[pi].s2 };

        ld[0] = &tls_dw.ldat[sp[0]];
        ld[1] = sp[1] == REFERENCE_SAMPLE ? &tls_dw.pseudo_sample : &tls_dw.ldat[sp[1]];
        
        tls_dw.metrics.total++;
        tls_dw.pair_stats[pi].total++;

        /* retrieve diff_state from pre-computed array */
        diff_state = it->cur->state;
        cache_hit = (it->cur->state != STATE_UNKNOWN);
        if (diff_state == STATE_UNKNOWN) {
            /* the pre-scanning didn't manage to get this from the
               cache. must be calculated here. */
            for (i = 0; i != 2; ++i) {
                if (! ld[i]->init.base_ct) {
                    ld[i]->base_ct = pileup_current_basecalls(sp[i]);
                    ld[i]->init.base_ct = 1;
                }
            }

            diff_state =
                dir_cache_calc_state(ld[0]->base_ct.ct_filt,
                                     ld[1]->base_ct.ct_filt,
                                     &ld[0]->dist,
                                     &ld[1]->dist);
        }
        
        tls_dw.pair_stats[pi].dist_count[diff_state]++;
        tls_dw.pair_stats[pi].cache_was_set += cache_hit;
        tls_dw.metrics.cache_was_set += cache_hit;

        if (diff_state == CHANGED) {

            struct pileup_locus_info pli;
            pileup_current_info(&pli);

            /* Finish sampling and do full distance marginal estimation */
            for (i = 0; i != 2; ++i) {
                if (! ld[i]->init.bqs_ct) {
                    pileup_current_bqs(sp[i], &ld[i]->bqs_ct, &ld[i]->n_bqs_ct);
                    ld[i]->init.bqs_ct = 1;
                }
                dir_points_update_alpha(ld[i]->base_ct.ct_filt, NULL, &ld[i]->dist);
                dir_points_fill(&ld[i]->dist);
                dir_weights_update_terms(ld[i]->bqs_ct, ld[i]->n_bqs_ct, &ld[i]->dist);
                if (! is_pseudo_sample(&ld[i]->dist))
                    dir_weights_fill(&ld[i]->dist);
            }
            compute_wsq_dist((const double *)ld[0]->dist.data, ld[0]->dist.weights,
                             (const double *)ld[1]->dist.data, ld[1]->dist.weights,
                             g_ld_par.max_sample_points,
                             tls_dw.square_dist_buf,
                             tls_dw.weights_buf);

            double test_quantile = 1.0 - g_ld_par.post_confidence;
            quantile_vals_t test_quantile_values;
            /* Compute the test distance quantile (relative to the
               dist, not squared dist) */
            compute_dist_wgt_quantiles(tls_dw.square_dist_buf,
                                       tls_dw.weights_buf,
                                       g_ld_par.max_sample_points,
                                       &test_quantile,
                                       1,
                                       test_quantile_values);

            if (test_quantile_values[0] > g_ld_par.min_dirichlet_dist) {
                ++tls_dw.pair_stats[pi].confirmed_changed;
                ld[0]->confirmed_changed = 1;
                ld[1]->confirmed_changed = 1;

                if (out_buf) {
                    compute_dist_wgt_quantiles(tls_dw.square_dist_buf,
                                               tls_dw.weights_buf,
                                               g_ld_par.max_sample_points,
                                               g_ld_par.quantiles,
                                               g_ld_par.n_quantiles,
                                               tls_dw.dist_quantile_values);            
                    
                    /* quantile values are in squared distance terms, and
                       in the [0, 1] x 4 space of points.  The maximum
                       distance between such points is sqrt(2.0).  We want
                       to re-scale it to be 1. */
                    unsigned q;
                    for (q = 0; q != g_ld_par.n_quantiles; ++q)
                        tls_dw.dist_quantile_values[q] *= ONE_OVER_SQRT2;
                    
                    struct pileup_locus_info pli;
                    pileup_current_info(&pli);
                    print_distance_quantiles(pli.refname,
                                             pli.pos + 1,
                                             pli.refbase,
                                             pi,
                                             tls_dw.dist_quantile_values, 
                                             out_buf);
                }
            }
        }
        
    }
}


void
print_indel_distance_quantiles(size_t pair_index,
                               double *dist_quantile_values,
                               struct indel_pair_count *events,
                               size_t n_events,
                               struct managed_buf *mb)
{
    unsigned 
        s1 = bam_sample_pairs.m[pair_index].s1,
        s2 = bam_sample_pairs.m[pair_index].s2;
    
    struct locus_data *ld[] = {
        &tls_dw.ldat[s1],
        s2 == REFERENCE_SAMPLE ? &tls_dw.pseudo_sample : &tls_dw.ldat[s2]
    };


    unsigned space = (2 * MAX_LABEL_LEN) + 3 + 100 + (10 * MAX_NUM_QUANTILES);
    ALLOC_GROW(mb->buf, mb->size + space, mb->alloc);

    struct pileup_locus_info pli;
    pileup_current_info(&pli);

    mb->size += sprintf(mb->buf + mb->size, 
                        "%s\t%s\t%s\t%c\t%u", 
                        bam_samples.m[s1].label,
                        s2 == REFERENCE_SAMPLE ? "REF" : bam_samples.m[s2].label,
                        pli.refname, 
                        pli.refbase, 
                        pli.pos + 1);
    
    /* indel distance quantiles */
    for (size_t q = 0; q != g_ld_par.n_quantiles; ++q)
        mb->size += sprintf(mb->buf + mb->size, "\t%.4f", dist_quantile_values[q]);

    unsigned indel_space = n_events * (10 + 10);
    ALLOC_GROW(mb->buf, mb->size + indel_space, mb->alloc);
    
    /* indel event counts for each sample.  The first event is the
       match event.  All remaining ones are counts of each type of
       indel. */
    struct indel_pair_count *eb, *ee = events + n_events;
    unsigned s;
    for (s = 0; s != 2; ++s) {
        mb->size += sprintf(mb->buf + mb->size, "\t%u",
                            ld[s]->base_ct.n_match_lo_q +
                            ld[s]->base_ct.n_match_hi_q);
        for (eb = events; eb != ee; ++eb)
            mb->size += sprintf(mb->buf + mb->size, ",%u", eb->count[s]);
    }

    /* retrieve each indel one by one. represent the union of indels
       in both samples as a comma-separated list.  insertions start
       with a '+', deletions with a '-'.  The match state is
       represented by '@'.  For example: @,+TCC,+TC,-G */
    mb->size += sprintf(mb->buf + mb->size, "\t@");
    static char ins[] = "-+";
    for (eb = events; eb != ee; ++eb) {
        struct indel_seq *isq = pileup_current_indel_seq(&eb->indel);
        unsigned grow = 2 + strlen(isq->seq);
        ALLOC_GROW(mb->buf, mb->size + grow, mb->alloc);
        mb->size += sprintf(mb->buf + mb->size, ",%c%s", 
                            ins[(unsigned)isq->is_ins], isq->seq);
        free(isq);
    }
    
    /* optional pileup fields.  */
    if (g_ld_par.do_print_pileup) {
        struct pileup_data
            *lp0 = &ld[0]->sample_data,
            *lp1 = &ld[1]->sample_data;
        
        unsigned extra_space = 
            lp0->calls.size + lp0->quals.size
            + lp1->calls.size + lp1->quals.size
            + 50;
        
        ALLOC_GROW(mb->buf, mb->size + extra_space, mb->alloc);
        char *out = mb->buf + mb->size;
        out += sprintf(out, "\t%u\t%u\t%u\t", 
                       lp0->n_match_hi_q, 
                       lp0->n_match_lo_q,
                       lp0->n_indel);

        strncpy(out, lp0->calls.buf, lp0->calls.size);
        out += lp0->calls.size;
        *out++ = '\t';

        strncpy(out, lp0->quals.buf, lp0->quals.size);
        out += lp0->quals.size;

        out += sprintf(out, "\t%u\t%u\t%u\t",
                       lp1->n_match_hi_q,
                       lp1->n_match_lo_q,
                       lp1->n_indel);

        strncpy(out, lp1->calls.buf, lp1->calls.size);
        out += lp1->calls.size;
        *out++ = '\t';

        strncpy(out, lp1->quals.buf, lp1->quals.size);
        out += lp1->quals.size;
        mb->size = out - mb->buf;
    }
    mb->buf[mb->size++] = '\n';
}


// print out all next distance quantiles for indels
void
indel_distance_quantiles_aux(struct managed_buf *buf)
{
    struct pileup_locus_info pli;
    pileup_current_info(&pli);

    struct indel_pair_count *indel_pairs = NULL;
    struct locus_data *ld[2];
    for (size_t pi = 0; pi != bam_sample_pairs.n; ++pi) {
        unsigned s[] = { bam_sample_pairs.m[pi].s1, bam_sample_pairs.m[pi].s2 };
        
        ld[0] = &tls_dw.ldat[s[0]];
        ld[1] = s[1] == REFERENCE_SAMPLE ? &tls_dw.pseudo_sample : &tls_dw.ldat[s[1]];
        
        unsigned i;
        for (i = 0; i != 2; ++i)
            if (! ld[i]->init.indel_ct) {
                pileup_current_indels(s[i], &ld[i]->indel_ct, &ld[i]->n_indel_ct);
                ld[i]->init.indel_ct = 1;
            }

        if (ld[0]->n_indel_ct == 0
            && ld[1]->n_indel_ct == 0
            && g_ld_par.min_dirichlet_dist > 0) continue;
        
        for (i = 0; i != 2; ++i) {
            if (! ld[i]->init.base_ct) {
                ld[i]->base_ct = pileup_current_basecalls(s[i]);
                ld[i]->init.base_ct = 1;
            }
            if (! ld[i]->init.sample_data) {
                pileup_current_data(s[i], &ld[i]->sample_data);
                ld[i]->init.sample_data = 1;
            }
        }


        
        /* traverse the sorted indel counts in tandem */
        unsigned n_indel_pairs = 1;
        pileup_make_indel_pairs(ld[0]->indel_ct, ld[0]->n_indel_ct,
                                ld[1]->indel_ct, ld[1]->n_indel_ct,
                                &indel_pairs, &n_indel_pairs);
        
        /* need at least some indels */
        if (n_indel_pairs == 0) continue;
        
        /* if # of CIGAR match reads > 0 for either sample, we need
           one extra 'event' */
        unsigned n_match[] = { 
            ld[0]->base_ct.n_match_lo_q + ld[0]->base_ct.n_match_hi_q,
            ld[1]->base_ct.n_match_lo_q + ld[1]->base_ct.n_match_hi_q
        };

        /* in order to compare the two dirichlet posteriors, we
           consider that the CIGAR 'match' event exists if either one
           of the bam_samples has non-zero counts.  If both are zero,
           though, we don't include the CIGAR match event in the
           dirichlet's at all. */
        unsigned has_match_reads = (n_match[0] != 0 || n_match[1] != 0);
        unsigned n_events = n_indel_pairs + has_match_reads;

        double *alpha[] = { malloc(n_events * sizeof(double)),
                            malloc(n_events * sizeof(double))
        };
        size_t buf_sz = n_events * g_ld_par.max_sample_points;
        double *points[] = { malloc(buf_sz * sizeof(double)),
                             malloc(buf_sz * sizeof(double))
        };

        unsigned c;
        double *p, *pe;
        for (i = 0; i != 2; ++i) {
            for (c = 0; c != n_indel_pairs; ++c)
                alpha[i][c] = indel_pairs[c].count[i] + g_ld_par.indel_prior_alpha;
            if (has_match_reads)
                alpha[i][n_indel_pairs] = n_match[i] + g_ld_par.indel_prior_alpha;

            pe = points[i] + buf_sz;
            for (p = points[i]; p != pe; p += n_events)
                gsl_ran_dirichlet(tls_dw.randgen, n_events, alpha[i], p);
        }

        
        compute_square_dist(points[0], points[1], 
                            g_ld_par.max_sample_points, n_events,
                            tls_dw.square_dist_buf);

        double test_quantile = 1.0 - g_ld_par.post_confidence, test_quantile_value;
            
        // compute distance quantiles
        compute_marginal_quantiles(tls_dw.square_dist_buf,
                                   g_ld_par.max_sample_points,
                                   1, /* one dimensional */
                                   0, /* use the first dimension */
                                   &test_quantile,
                                   1, /* evaluate only one quantile */
                                   &test_quantile_value);

        test_quantile_value = sqrt(test_quantile_value);
        if (test_quantile_value >= g_ld_par.min_dirichlet_dist) {
            compute_marginal_quantiles(tls_dw.square_dist_buf,
                                       g_ld_par.max_sample_points,
                                       1, /* one dimensional */
                                       0, /* use the first dimension */
                                       g_ld_par.quantiles,
                                       g_ld_par.n_quantiles,
                                       tls_dw.dist_quantile_values);
            unsigned q;
            for (q = 0; q != g_ld_par.n_quantiles; ++q)
                tls_dw.dist_quantile_values[q] = 
                    sqrt(tls_dw.dist_quantile_values[q]) * ONE_OVER_SQRT2;

            print_indel_distance_quantiles(pi, 
                                           tls_dw.dist_quantile_values,
                                           indel_pairs, n_indel_pairs, buf);
        }
        free(points[0]);
        free(points[1]);
        free(alpha[0]);
        free(alpha[1]);
    }
    free(indel_pairs);
}


/* compute and print base composition marginal quantiles. */
void
comp_quantiles_aux(struct managed_buf *comp_buf)
{
    comp_quantile_vals_t comp_quantile_values;
    double comp_means[NUM_NUCS];

    struct pileup_locus_info ploc;
    ploc.pos = UINT_MAX; /* signal that this isn't initialized yet */

    unsigned d, s;
    for (s = 0; s != bam_samples.n; ++s) {
        struct locus_data *ld = &tls_dw.ldat[s];
        
        if (ld->confirmed_changed) {

            dir_points_fill(&ld->dist);
            de_permute_points(&ld->dist);
            dir_weights_fill(&ld->dist);
            
            compute_marginal_wgt_quantiles((double *)ld->dist.data,
                                           ld->dist.weights,
                                           ld->dist.n_points,
                                           NUM_NUCS,
                                           g_ld_par.quantiles,
                                           g_ld_par.n_quantiles,
                                           (quantile_vals_t *)&comp_quantile_values[0][0]);
            for (d = 0; d != NUM_NUCS; ++d)
                comp_means[d] = 
                    compute_marginal_mean((double *)ld->dist.data,
                                          ld->dist.weights,
                                          ld->dist.n_points,
                                          NUM_NUCS,
                                          d);

            if (! ld->init.sample_data) {
                pileup_current_data(s, &ld->sample_data);
                ld->init.sample_data = 1;
            }
            
            /* just-in-time initialization of ploc */
            if (ploc.pos == UINT_MAX)
                pileup_current_info(&ploc);
            
            print_basecomp_quantiles(comp_quantile_values,
                                     comp_means,
                                     g_ld_par.n_quantiles,
                                     bam_samples.m[s].label,
                                     &ploc,
                                     &ld->sample_data,
                                     comp_buf);
        }
    }
}


/* order by pk */
static int
pair_key_type_less(struct pos_pair_key a,
                   struct pos_pair_key b)
{
    return a.pk.key < b.pk.key
        || (a.pk.key == b.pk.key
            && a.pk.is_ref_change < b.pk.is_ref_change);
}

KSORT_INIT(key_type_sort, struct pos_pair_key, pair_key_type_less);


static int
pair_pos_less(struct pos_pair_key a,
              struct pos_pair_key b)
{
    return cmp_contig_pos(a.pos, b.pos) == -1;
}


KSORT_INIT(pos_sort, struct pos_pair_key, pair_pos_less);


/* prescan all input to initialize keys.  then, sort by key, allowing
   a single computation per distinct key value in sequence.  then,
   populate 'state' variable.  finally, sort by position again. */
static struct pos_pair_iter *
prescan_input()
{
    struct pos_pair_iter *dist_class =
        malloc(bam_sample_pairs.n * sizeof(struct pos_pair_iter));
    
    unsigned s, pi;
    for (pi = 0; pi != bam_sample_pairs.n; ++pi)
        dist_class[pi] = (struct pos_pair_iter){ NULL, NULL, 0, 0 };

    /* pre-scan all input */
    struct locus_data *ld[2];
    while (pileup_next_pos()) {
        for (s = 0; s != bam_samples.n; ++s) {
            reset_locus_data(&tls_dw.ldat[s]);
            tls_dw.ldat[s].base_ct = pileup_current_basecalls(s);
        }
        reset_locus_data(&tls_dw.pseudo_sample);
        tls_dw.pseudo_sample.base_ct =
            pileup_current_basecalls(PSEUDO_SAMPLE);

        struct contig_pos pos = pileup_current_pos();
        for (pi = 0; pi != bam_sample_pairs.n; ++pi) {
            unsigned sp[] = { bam_sample_pairs.m[pi].s1,
                              bam_sample_pairs.m[pi].s2 };
            ld[0] = &tls_dw.ldat[sp[0]];
            ld[1] = sp[1] == REFERENCE_SAMPLE
                ? &tls_dw.pseudo_sample :
                &tls_dw.ldat[sp[1]];
            struct pos_pair_iter *pp = &dist_class[pi];
            ALLOC_GROW(pp->buf, pp->size + 1, pp->alloc);
            pp->buf[pp->size].pos = pos;
            pp->buf[pp->size].pk =
                dir_cache_classify_alphas(ld[0]->base_ct.ct_filt,
                                          ld[1]->base_ct.ct_filt);
            ++pp->size;
        }
    }

    /* sort by key value and key type. populate 'state' variable from
       cache if possible. */
    unsigned i;
    struct pos_pair_key pps = { .pk = { .is_valid = 0 } }; /* previous pk */
    for (pi = 0; pi != bam_sample_pairs.n; ++pi) {
        struct pos_pair_iter *pp = &dist_class[pi];
        ks_introsort(key_type_sort, pp->size, pp->buf);
        for (i = 0; i != pp->size; ++i) {
            struct pos_pair_key *ps = &pp->buf[i];
            if (ps->pk.is_valid) {
                if (pps.pk.is_valid && pps.pk.key == ps->pk.key)
                    ps->state = pps.state;
                else
                    ps->state = dir_cache_try_get_state(ps->pk);
            } else
                ps->state = STATE_UNKNOWN;
            pps = *ps;
        }
    }

    /* re-sort back to position based sorting. */
    for (pi = 0; pi != bam_sample_pairs.n; ++pi) {
        struct pos_pair_iter *pp = &dist_class[pi];
        ks_introsort(pos_sort, pp->size, pp->buf);
        pp->cur = pp->buf;
    }
    pileup_reset_pos();
    
    return dist_class;
}



void
locus_diff_worker(const struct managed_buf *in_bufs,
                  unsigned more_input,
                  void *vsi,
                  struct managed_buf *out_bufs)
{
    struct timespec worker_start_time;
    clock_gettime(CLOCK_REALTIME, &worker_start_time);
    tls_dw.metrics.total = 0;
    tls_dw.metrics.cacheable = 0;
    tls_dw.metrics.cache_was_set = 0;

    unsigned i = 0;
    struct managed_buf 
        *dist_buf = g_ld_par.do_dist ? &out_bufs[i++] : NULL,
        *comp_buf = g_ld_par.do_comp ? &out_bufs[i++] : NULL,
        *indel_buf = g_ld_par.do_indel ? &out_bufs[i++] : NULL;

    struct managed_buf bam = { NULL, 0, 0 };
    unsigned s;
    struct bam_scanner_info *bsi = vsi;

    pileup_load_refseq_ranges(bsi);
    for (s = 0; s != bam_samples.n; ++s) {
        struct bam_stats *bs = &bsi->m[s];
        bam_inflate(&in_bufs[s], bs->chunks, bs->n_chunks, &bam);
        pileup_tally_stats(bam, bsi, s);
    }
    if (bam.buf != NULL) free(bam.buf);

    if (! more_input)
        pileup_final_input();

    /* we need all three types of data, so must prepare all of it. */
    for (s = 0; s != bam_samples.n; ++s) {
        pileup_prepare_bqs(s);
        pileup_prepare_basecalls(s);
        pileup_prepare_indels(s);
    }

    /* zero out the pair_stats */
    unsigned pi;
    for (pi = 0; pi != bam_sample_pairs.n; ++pi)
        memset(&tls_dw.pair_stats[pi], 0, sizeof(tls_dw.pair_stats[0]));

    /* allocate pre-computed distances. */
    struct pos_pair_iter *dist_class = prescan_input();

    while (pileup_next_pos()) {
        reset_locus_data(&tls_dw.pseudo_sample);
        for (s = 0; s != bam_samples.n; ++s)
            reset_locus_data(&tls_dw.ldat[s]);

        if (dist_buf || comp_buf)
            distance_quantiles_aux(dist_buf, dist_class);
        
        if (indel_buf)
            indel_distance_quantiles_aux(indel_buf);

        if (comp_buf)
            comp_quantiles_aux(comp_buf);

        for (pi = 0; pi != bam_sample_pairs.n; ++pi)
            dist_class[pi].cur++;
    }   

    /* frees statistics that have already been used in one of the
       distance calculations. */
    pileup_clear_stats();
    accumulate_pair_stats(tls_dw.pair_stats);

    for (pi = 0; pi != bam_sample_pairs.n; ++pi)
        free(dist_class[pi].buf);
    free(dist_class);
}
 
    
void
locus_diff_offload(void *par, const struct managed_buf *bufs)
{
    struct locus_diff_offload_par *ol = (struct locus_diff_offload_par *)par;
    unsigned i = 0;
    if (ol->dist_fh) {
        fwrite(bufs[i].buf, 1, bufs[i].size, ol->dist_fh), i++;
        fflush(ol->dist_fh);
    }

    if (ol->comp_fh) {
        fwrite(bufs[i].buf, 1, bufs[i].size, ol->comp_fh), i++;
        fflush(ol->comp_fh);
    }

    if (ol->indel_fh) {
        fwrite(bufs[i].buf, 1, bufs[i].size, ol->indel_fh), i++;
        fflush(ol->indel_fh);
    }
}
