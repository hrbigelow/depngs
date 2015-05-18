#define __STDC_LIMIT_MACROS
#include <stdint.h>

#include "dist_worker.h"
#include "sampling.h"
#include "yepLibrary.h"

#include <gsl/gsl_math.h>
#include <gsl/gsl_randist.h>
#include <time.h>
#include <pthread.h>


extern "C" {
#include "locus.h"
#include "dirichlet_points_gen.h"
#include "dirichlet_diff_cache.h"
#include "khash.h"
#include "range_line_reader.h"
}

#define MIN(a, b) ((a) < (b) ? (a) : (b))

struct pair_dist_stats {
    size_t total; /* total number of loci pairs for which at least one
                     sample has coverage */
    size_t dist_count[5]; /* number of loci in 'total' that fall into
                             each of the categories */
    size_t confirmed_changed; /* number of dist_count[CHANGED] loci
                                 that are further confirmed by
                                 weighted sampling. */
    size_t cacheable; /* number of loci that are primary or secondary
                         cacheable, based on the 8 alpha values.  */
    size_t cache_was_set; /* number of cacheable loci in which the
                             cache was already set. */
};

static timespec start_time;
static double posterior_confidence;
static double min_dirichlet_dist;
static unsigned max_sample_points;

#define PSEUDO_SAMPLE -1

/* describes all pairs of samples */
static struct {
    struct { 
        int s1, s2; 
        struct pair_dist_stats stats;
    } *p;
    unsigned n;
} sample_pairs;

static pthread_mutex_t pair_stats_mtx = PTHREAD_MUTEX_INITIALIZER;


#define MAX_LABEL_LEN 100
/* attributes intrinsic to one sample */
static struct {
    struct {
        char *file;
        char label[MAX_LABEL_LEN + 1];
        struct nucleotide_stats nuc_stats;
        FILE *fh;
    } *atts;
    unsigned n;
} samples;

static struct {
    unsigned do_print_pileup;
    unsigned do_dist, do_comp, do_indel;
} worker_options;


static struct {
    struct range_line_reader_par *reader_buf;
    void **reader_par;
    unsigned n_readers;
    struct pair_ordering_range *ranges;
    struct pair_ordering global_read_start;
    unsigned n_ranges;
    struct dist_worker_input *w;
    void **wp;
    unsigned n_threads;
    struct dist_worker_offload_par offload_par;
} thread_params;

static double quantiles[MAX_NUM_QUANTILES];
static unsigned n_quantiles;


KHASH_MAP_INIT_STR(remap, int)

void init_sample_attributes(const char *samples_file,
                            const char *fastq_type,
                            khash_t(remap) *map);

void init_sample_pairs(const char *pair_file, khash_t(remap) *map);
void init_quantiles(const char *csv, double *vals, unsigned *n_vals);


void alloc_pileup_locus(struct locus_sampling *ls)
{
    ls->locus = PileupSummary();
    ls->is_next = 0;
    ls->confirmed_changed = 0;
    alloc_distrib_points(&ls->distp);
    ls->locus_ord = (struct pair_ordering){ 0, 0 };
}


void free_pileup_locus(struct locus_sampling *ls)
{
    free_distrib_points(&ls->distp);
}



void dist_worker_init(double _post_confidence, 
                      double _min_dirichlet_dist,
                      unsigned _max_sample_points,
                      const char *samples_file,
                      const char *sample_pairs_file,
                      const char *fastq_type,
                      const char *quantiles_string,
                      unsigned do_dist,
                      unsigned do_comp,
                      unsigned do_indel,
                      unsigned do_print_pileup)
{
    clock_gettime(CLOCK_REALTIME, &start_time);
    if (_post_confidence < 0.8 || _post_confidence > 0.999999)
    {
        fprintf(stderr,
                "Error: dist_worker_init: posterior confidence of %g is a bad value\n",
                _post_confidence);
        exit(1);
    }
    posterior_confidence = _post_confidence;
    if (_min_dirichlet_dist <= 0 || _min_dirichlet_dist >= 1)
    {
        fprintf(stderr,
                "Error: dist_worker_init: min_dirichlet_dist of %g is a bad value.\n"
                "Should be in [0, 1]\n", _min_dirichlet_dist);
        exit(1);
    }
    min_dirichlet_dist = _min_dirichlet_dist;
    max_sample_points = _max_sample_points;

    khash_t(remap) *map = kh_init(remap);
    init_sample_pairs(sample_pairs_file, map);
    init_sample_attributes(samples_file, fastq_type, map);

    khiter_t k;
    for (k = kh_begin(map); k != kh_end(map); ++k)
        if (kh_exist(map, k)) free((char *)kh_key(map, k));
    kh_destroy(remap, map);

    init_quantiles(quantiles_string, quantiles, &n_quantiles);

    worker_options.do_print_pileup = do_print_pileup;
    worker_options.do_dist = do_dist;
    worker_options.do_comp = do_comp;
    worker_options.do_indel = do_indel;
}


void dist_worker_free()
{
    free(sample_pairs.p);
    unsigned s;
    for (s = 0; s != samples.n; ++s)
    {
        free(samples.atts[s].file);
        fclose(samples.atts[s].fh);
    }
    free(samples.atts);
}



void init_sample_attributes(const char *samples_file,
                            const char *fastq_type,
                            khash_t(remap) *map)
{
    char jpd[1000], label[1000], pileup[1000];
    samples.n = kh_size(map);
    unsigned s, alloc = 0;
    ALLOC_GROW_TYPED(samples.atts, samples.n, alloc);
    khiter_t k;
    FILE *samples_fh = open_if_present(samples_file, "r");
    while (! feof(samples_fh))
    {
        (void)fscanf(samples_fh, "%s\t%s\t%s\n", label, jpd, pileup);
        if ((k = kh_get(remap, map, label)) != kh_end(map))
        {
            s = kh_val(map, k);
            samples.atts[s].file = strdup(pileup);
            nucleotide_stats_initialize(jpd, &samples.atts[s].nuc_stats);
            if (strlen(label) + 1 > sizeof(samples.atts[s].label)) 
            {
                fprintf(stderr, "%s: error: sample label string must be "
                        "less than %Zu characters\n", __func__,
                        sizeof(samples.atts[s].label));
                exit(1);
            }
            strcpy(samples.atts[s].label, label);
            samples.atts[s].fh = fopen(pileup, "r");
            if (! samples.atts[s].fh)
            {
                fprintf(stderr, "%s: error: couldn't open pileup input file %s\n",
                        __func__, pileup);
                exit(1);
            }
        }
    }
    fclose(samples_fh);

    /* initialize fastq type based on the first */
    int offset;
    if (fastq_type) offset = fastq_type_to_offset(fastq_type);
    else
    {
        size_t sz = 1024 * 1024 * 64;
        char *buf = (char *)malloc(sz);

        /* using just the first sample file here. this assumes all
           pileup files mentioned have the same offset. */
        offset = fastq_offset(samples.atts[0].file, buf, sz);
        free(buf);
    }

    if (offset == -1)
    {
        fprintf(stderr, "Could not determine fastq type of this pileup file.\n");
        exit(1);
    }

    PileupSummary::set_offset(offset);
}


/* parse a sample pairs file. */
void init_sample_pairs(const char *pair_file, khash_t(remap) *map)
{
    unsigned idx = 0;
    int ret;
    FILE *sample_pair_fh = open_if_present(pair_file, "r");
    char label[2][MAX_LABEL_LEN];
    khiter_t k1, k2;
    unsigned alloc = 0;

    /* initialize the keyword PSEUDO as the reference sample */
    char pseudo_key[] = "PSEUDO";
    k1 = kh_put(remap, map, pseudo_key, &ret);
    assert(ret == 0 || ret == 1);
    kh_val(map, k1) = -1;

    unsigned p = 0;
    while (! feof(sample_pair_fh))
    {
        (void)fscanf(sample_pair_fh, "%s\t%s\n", label[0], label[1]);
        if ((k1 = kh_get(remap, map, label[0])) == kh_end(map))
        {
            k1 = kh_put(remap, map, strdup(label[0]), &ret);
            assert(ret == 0 || ret == 1);
            kh_val(map, k1) = idx++;
        }
        if ((k2 = kh_get(remap, map, label[1])) == kh_end(map))
        {
            k2 = kh_put(remap, map, strdup(label[1]), &ret);
            assert(ret == 0 || ret == 1);
            kh_val(map, k2) = idx++;
        }
        ALLOC_GROW_TYPED(sample_pairs.p, p + 1, alloc);
        sample_pairs.p[p].s1 = kh_val(map, k1);
        sample_pairs.p[p].s2 = kh_val(map, k2);
        ++p;
    }
    sample_pairs.n = p;
    fclose(sample_pair_fh);

    k1 = kh_get(remap, map, pseudo_key);
    kh_del(remap, map, k1);
}


/* parse a comma-separated list of values into an array. check
   that the values are between 0 and 1 and increasing. */
void init_quantiles(const char *csv, double *vals, unsigned *n_vals)
{
    *n_vals = 0;
    errno = 0;
    double pv = -1;
    unsigned err = 0;
    char *loc;
    while (*csv != '\0' && *n_vals != MAX_NUM_QUANTILES)
    {    
        if (*csv == ',') ++csv;
        *vals = strtod(csv, &loc);
        csv = loc;
        err = (err || *vals < 0 || *vals > 1 || *vals < pv);
        pv = *vals;
        ++vals;
        ++(*n_vals);
    }
    err = (err || *csv != '\0' || *n_vals == 0);
    if (errno || err)
    {
        fprintf(stderr, "%s: The input is supposed to be a comma-separated"
                "list of between 1 and %u increasing numbers in the [0, 1] range\n"
                "Input was: %s\n",
                __func__, MAX_NUM_QUANTILES, csv);
        exit(1);
    }
}


#define N_STATS_CATEGORIES \
    sizeof(sample_pairs.p[0].stats.dist_count) \
    / sizeof(sample_pairs.p[0].stats.dist_count[0])

void accumulate_pair_stats(struct pair_dist_stats *stats)
{
    pthread_mutex_lock(&pair_stats_mtx);
    unsigned s, p;
    for (p = 0; p != sample_pairs.n; ++p)
    {
        sample_pairs.p[p].stats.total += stats[p].total;
        sample_pairs.p[p].stats.confirmed_changed += stats[p].confirmed_changed;
        sample_pairs.p[p].stats.cacheable += stats[p].cacheable;
        sample_pairs.p[p].stats.cache_was_set += stats[p].cache_was_set;

        for (s = 0; s != N_STATS_CATEGORIES; ++s)
            sample_pairs.p[p].stats.dist_count[s] += stats[p].dist_count[s];
    }

    pthread_mutex_unlock(&pair_stats_mtx);
}

#define NUM_EXTRA_BUFS_PER_THREAD 500

struct thread_queue *dist_worker_tq_init(const char *query_range_file,
                                         unsigned n_threads,
                                         unsigned n_readers,
                                         unsigned long max_input_mem,
                                         FILE *dist_fh,
                                         FILE *comp_fh,
                                         FILE *indel_fh)
{
    thread_params.n_threads = n_threads;
    thread_params.w = (struct dist_worker_input *)
        malloc(n_threads * sizeof(struct dist_worker_input));
    
    thread_params.wp = (void **)malloc(n_threads * sizeof(void *));
    struct locus_sampling *ls;
    unsigned s, t;
    for (t = 0; t != n_threads; ++t)
    {
        thread_params.w[t].randgen = gsl_rng_alloc(gsl_rng_taus);
        alloc_pileup_locus(&thread_params.w[t].pseudo_sample);
        thread_params.w[t].lslist = new struct locus_sampling[samples.n];
        for (s = 0; s != samples.n; ++s)
        {
            ls = &thread_params.w[t].lslist[s];
            alloc_pileup_locus(ls);
        }            
        thread_params.w[t].pair_stats =
            (struct pair_dist_stats *)calloc(sample_pairs.n, sizeof(struct pair_dist_stats));
        thread_params.w[t].square_dist_buf = (double *)malloc(sizeof(double) * max_sample_points);
        thread_params.w[t].weights_buf = (double *)malloc(sizeof(double) * max_sample_points);
        thread_params.w[t].do_print_progress = (t == n_threads - 1); /* last thread prints progress */
        thread_params.w[t].bep.points_hash_frozen = 0;
        thread_params.wp[t] = &thread_params.w[t];
    }

    size_t scan_thresh_size = 1e5;
    file_bsearch_init(init_locus, scan_thresh_size);

    thread_params.reader_buf = 
        (struct range_line_reader_par *)malloc(n_readers * sizeof(struct range_line_reader_par));

    thread_params.reader_par = (void **)malloc(n_readers * sizeof(void *));

    unsigned r;
    for (r = 0; r != n_readers; ++r)
    if (query_range_file)
        thread_params.ranges = 
            parse_query_ranges(query_range_file, &thread_params.n_ranges);
    else
    {
        /* simply set the 'query' to the largest span possible */
        thread_params.n_ranges = 1;
        thread_params.ranges = (struct pair_ordering_range *)
            malloc(sizeof(struct pair_ordering_range));
        thread_params.ranges[0] = 
            (struct pair_ordering_range){ { 0, 0 }, { SIZE_MAX, SIZE_MAX - 1 } };
    }

    for (r = 0; r != n_readers; ++r)
    {
        thread_params.reader_buf[r] = (struct range_line_reader_par){
            (struct file_bsearch_index *)malloc(samples.n 
                                                * sizeof(struct file_bsearch_index)),
            samples.n,
            thread_params.ranges, 
            thread_params.ranges + thread_params.n_ranges,
            { 0, 0 },
            { 0, 0 },
            init_locus,
            1
        };
        for (s = 0; s != samples.n; ++s)
            thread_params.reader_buf[r].ix[s] =
                file_bsearch_make_index(samples.atts[s].file);

        thread_params.reader_par[r] = &thread_params.reader_buf[r];
    }
    
    thread_params.offload_par = { dist_fh, comp_fh, indel_fh };

    enum YepStatus status = yepLibrary_Init();
    assert(status == YepStatusOk);

    /* To avoid a stall, n_extra / n_threads should be greater than
       Max(work chunk time) / Avg(work chunk time). */
    size_t n_extra = n_threads * NUM_EXTRA_BUFS_PER_THREAD;
    size_t n_output_files = 
        (dist_fh ? 1 : 0) + (comp_fh ? 1 : 0) + (indel_fh ? 1 : 0);

    thread_queue_reader_t reader = { rl_reader, rl_scanner, rl_get_start, rl_set_start };
    thread_params.global_read_start = (struct pair_ordering){ 0, 0 };

    struct thread_queue *tqueue =
        thread_queue_init(reader, thread_params.reader_par,
                          dist_worker, thread_params.wp,
                          dist_offload, &thread_params.offload_par,
                          &thread_params.global_read_start,
                          n_threads, n_extra, n_readers, samples.n,
                          n_output_files, max_input_mem);

    return tqueue;
}


void dist_worker_tq_free()
{
    unsigned r, s;
    for (r = 0; r != thread_params.n_readers; ++r)
    {
        for (s = 0; samples.n; ++s)
            file_bsearch_index_free(thread_params.reader_buf[r].ix[s]);

        free(thread_params.reader_buf[r].ix);
    }

    free(thread_params.reader_buf);
    free(thread_params.reader_par);
    free(thread_params.ranges);

    unsigned t;
    for (t = 0; t != thread_params.n_threads; ++t)
    {
        gsl_rng_free(thread_params.w[t].randgen);
        free_pileup_locus(&thread_params.w[t].pseudo_sample);
        for (s = 0; s != samples.n; ++s)
            free_pileup_locus(&thread_params.w[t].lslist[s]);

        delete[] thread_params.w[t].lslist;
        free(thread_params.w[t].pair_stats);
        free(thread_params.w[t].square_dist_buf);
        free(thread_params.w[t].weights_buf);
    }
    free(thread_params.w);
    free(thread_params.wp);
}

void print_pair_stats(const char *stats_file)
{
    FILE *fh = open_if_present(stats_file, "w");
    pthread_mutex_lock(&pair_stats_mtx);
    fprintf(fh, 
            "%s\t%s\t%s\t%s\t%s\t%s", 
            "sample1", "sample2", "total", "cacheable",
            "cache_was_set", "confirmed_changed");

    unsigned s;
    for (s = 0; s != N_STATS_CATEGORIES; ++s)
        fprintf(fh, "\t%s", fuzzy_state_strings[s]);
    fprintf(fh, "\n");

    unsigned p;
    for (p = 0; p != sample_pairs.n; ++p)
    {
        /* Print out statistics */
        fprintf(fh, "%s\t%s", 
                samples.atts[sample_pairs.p[p].s1].label,
                samples.atts[sample_pairs.p[p].s2].label);
        fprintf(fh, "\t%zu", sample_pairs.p[p].stats.total);
        fprintf(fh, "\t%zu", sample_pairs.p[p].stats.cacheable);
        fprintf(fh, "\t%zu", sample_pairs.p[p].stats.cache_was_set);
        fprintf(fh, "\t%zu", sample_pairs.p[p].stats.confirmed_changed);

        for (s = 0; s != N_STATS_CATEGORIES; ++s)
            fprintf(fh, "\t%zu", sample_pairs.p[p].stats.dist_count[s]);
        fprintf(fh, "\n");
    }
    fclose(fh);
    pthread_mutex_unlock(&pair_stats_mtx);
}

typedef double COMP_QV[NUM_NUCS][MAX_NUM_QUANTILES];

void print_basecomp_quantiles(COMP_QV quantile_values,
                              const double *means,
                              size_t n_quantiles,
                              const char *label_string,
                              struct locus_sampling *ls, 
                              struct managed_buf *mb)
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

    unsigned grow = NUM_NUCS * (sizeof(line_label) + 30 + (11 * MAX_NUM_QUANTILES));
    ALLOC_GROW_TYPED(mb->buf, mb->size + grow, mb->alloc);

    static const char dimension_labels[] = "ACGT";
    for (d = 0; d != NUM_NUCS; ++d)
    {
        mb->size += sprintf(mb->buf + mb->size,
                            "%s\t%c\t%Zu\t%10.8f", line_label, 
                            dimension_labels[d], mean_rank_order[d], means[d]);
        
        for (q = 0; q != n_quantiles; ++q)
            mb->size += sprintf(mb->buf + mb->size, "\t%10.8f", quantile_values[d][q]);

        mb->size += sprintf(mb->buf + mb->size, "\n");
    }
}


/* populates square_dist_buf with squares of euclidean distances
   between points1 and points2 in the normalized [0,1] x 4 space.
   here, the maximum squared distance will be 2 (e.g. (1,0,0,0) and
   (0,1,0,0).  This will yield a maximum distance of sqrt(2).
   populates weights_buf with product of weights1 and weights2. */
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
                             const double *quantiles,
                             size_t n_quantiles,
                             double *dist_quantile_values)
{
    compute_marginal_wquantiles(square_dist_buf, weights_buf, n_points, 1, 0,
                                quantiles, n_quantiles,
                                dist_quantile_values);

    unsigned q;
    for (q = 0; q != n_quantiles; ++q) 
        dist_quantile_values[q] = sqrt(dist_quantile_values[q]);
}


/* print out distance quantiles. use a pseudo-sample for the second
   one if its pair index is equal to PSEUDO_SAMPLE */
void print_distance_quantiles(const char *contig,
                              size_t position,
                              struct dist_worker_input *dw,
                              size_t pair_index,
                              double *dist_quantile_values,
                              struct managed_buf *buf)
{
    int s1 = sample_pairs.p[pair_index].s1,
        s2 = sample_pairs.p[pair_index].s2;

    unsigned space = (2 * MAX_LABEL_LEN) + 3 + 100 + (10 * MAX_NUM_QUANTILES);

    ALLOC_GROW_TYPED(buf->buf, buf->size + space, buf->alloc);

    buf->size += sprintf(buf->buf + buf->size, 
                         "%s\t%s\t%s\t%Zu", 
                         samples.atts[s1].label,
                         s2 == PSEUDO_SAMPLE ? "REF" : samples.atts[s2].label,
                         contig, position);
    
    for (size_t q = 0; q != n_quantiles; ++q)
        buf->size += sprintf(buf->buf + buf->size, "\t%7.4f", dist_quantile_values[q]);

    locus_sampling *ls1 = &dw->lslist[s1], 
        *ls2 = s2 == PSEUDO_SAMPLE ? &dw->pseudo_sample : &dw->lslist[s2];

    if (worker_options.do_print_pileup)
    {
        unsigned extra_space = 
            ls1->locus.bases_raw.size
            + ls1->locus.quality_codes.size
            + ls2->locus.bases_raw.size
            + ls2->locus.quality_codes.size
            + 50;

        ALLOC_GROW_TYPED(buf->buf, buf->size + extra_space, buf->alloc);
        
        buf->size += sprintf(buf->buf + buf->size, 
                             "\t%Zu\t%s\t%s\t%Zu\t%s\t%s",
                             ls1->locus.read_depth,
                             ls1->locus.bases_raw.buf,
                             ls1->locus.quality_codes.buf,
                             ls2->locus.read_depth,
                             ls2->locus.bases_raw.buf,
                             ls2->locus.quality_codes.buf);
    }
    buf->size += sprintf(buf->buf + buf->size, "\n");
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



/* update 'ls' fields to be consistent with sd->current */
void update_pileup_locus(const struct nucleotide_stats *stats,
                         locus_sampling *ls)
{
    ls->locus.load_line(ls->current);
    ls->locus_ord = init_locus(ls->current);
    ls->locus.parse();
    nucleotide_stats_pack(stats, &ls->locus.counts);
    ls->distp.points.size = 0;
    ls->distp.weights.size = 0;
}


void init_pseudo_locus(struct locus_sampling *ls)
{
    ls->is_next = 1;
    ls->locus.read_depth = PSEUDO_DEPTH;
    ls->locus.read_depth_match = PSEUDO_DEPTH;
    ls->locus.read_depth_high_qual = PSEUDO_DEPTH;

    ls->distp.points.size = 0;
    ls->distp.weights.size = 0;

    ((struct points_gen_par *)ls->distp.pgen.points_gen_par)->post_counts = 
        &ls->locus.counts;

    ls->locus.bases_raw.size = 4;
    ALLOC_GROW_TYPED(ls->locus.bases_raw.buf, 
                     ls->locus.bases_raw.size,
                     ls->locus.bases_raw.alloc);
    strcpy(ls->locus.bases_raw.buf, "REF");

    ls->locus.quality_codes.size = 4;
    ALLOC_GROW_TYPED(ls->locus.quality_codes.buf, 
                     ls->locus.quality_codes.size,
                     ls->locus.quality_codes.alloc);
    strcpy(ls->locus.quality_codes.buf, "REF");

}

/* update ls->locus.counts with ultra-high depth of perfect 'nuc'
   basecalls. */
void update_pseudo_locus(char nuc, struct locus_sampling *ls)
{
    unsigned inuc = (unsigned)base_to_index(nuc);
    if (inuc == 4)
        ls->locus.counts.num_data = 0;
    else
    {
        ls->locus.counts.num_data = 1;
        ls->locus.counts.stats_index[0] = 
            encode_nucleotide(nuc, NUC_HIGHEST_QUALITY - 1, 0);
        unsigned b;
        for (b = 0; b != NUM_NUCS; ++b)
        {
            ls->locus.counts.stats[0].cpd[b] = b == inuc ? 1 : 0;
            ls->locus.base_counts_high_qual[b] = b == inuc ? PSEUDO_DEPTH : 0;
        }

        ls->locus.counts.stats[0].ct = PSEUDO_DEPTH;
    }
}

/* advance ls->current, then initialize if there is another locus */
void refresh_locus(const struct nucleotide_stats *stats,
                   locus_sampling *ls)
{
    assert(ls->current != ls->end);
    if ((ls->current = strchr(ls->current, '\n') + 1) == ls->end) 
        ls->is_next = 0;

    else
    {
        update_pileup_locus(stats, ls);
        ls->confirmed_changed = 0;
    }
}


void print_indel_distance_quantiles(const char *contig,
                                    size_t position,
                                    size_t pair_index,
                                    double *dist_quantile_values,
                                    indel_event *events,
                                    size_t n_events,
                                    locus_sampling *sd,
                                    struct managed_buf *mb)
{
    unsigned s1 = sample_pairs.p[pair_index].s1,
        s2 = sample_pairs.p[pair_index].s2;

    unsigned space = (2 * MAX_LABEL_LEN) + 3 + 100 + (10 * MAX_NUM_QUANTILES);
    ALLOC_GROW_TYPED(mb->buf, mb->size + space, mb->alloc);

    mb->size += sprintf(mb->buf + mb->size, 
                        "%s\t%s\t%s\t%Zu", 
                        samples.atts[s1].label,
                        samples.atts[s2].label,
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

    if (sd)
    {
        unsigned extra_space = 
            sd[s1].locus.bases_raw.size
            + sd[s1].locus.quality_codes.size
            + sd[s2].locus.bases_raw.size
            + sd[s2].locus.quality_codes.size
            + 50;

        ALLOC_GROW_TYPED(mb->buf, mb->size + extra_space, mb->alloc);
        mb->size += sprintf(mb->buf + mb->size,
                            "\t%Zu\t%s\t%s\t%Zu\t%s\t%s",
                            sd[s1].locus.read_depth,
                            sd[s1].locus.bases_raw.buf,
                            sd[s1].locus.quality_codes.buf,
                            sd[s2].locus.read_depth,
                            sd[s2].locus.bases_raw.buf,
                            sd[s2].locus.quality_codes.buf);
    }
    mb->size += sprintf(mb->buf + mb->size, "\n");
}

#define ONE_OVER_SQRT2 0.70710678118654752440

/* for each sample pair, calculate whether the loci differ above the
   given level and confidence. if a pair differs, set the
   confirmed_changed flag for each sample in the pair. if out_buf is
   not NULL, print out distance quantiles.  also generates sample
   points for each sample as needed, both for the preliminary test and
   more points for the final test */
void next_distance_quantiles_aux(struct dist_worker_input *dw,
                                 size_t gs,
                                 struct managed_buf *out_buf)
{
    locus_sampling *ls1, *ls2;
    size_t pi, i;
    unsigned counts[2][NUM_NUCS];
    enum fuzzy_state diff_state = AMBIGUOUS;

    unsigned cacheable, cache_was_set;

    for (pi = 0; pi != sample_pairs.n; ++pi)
    {
        int s1 = sample_pairs.p[pi].s1,
            s2 = sample_pairs.p[pi].s2;

        ls1 = &dw->lslist[s1];
        ls2 = s2 == PSEUDO_SAMPLE ? &dw->pseudo_sample : &dw->lslist[s2];
        
        dw->bep.dist[0] = &ls1->distp;
        dw->bep.dist[1] = &ls2->distp;

        for (i = 0; i != NUM_NUCS; ++i)
        {
            counts[0][i] = ls1->locus.base_counts_high_qual[i];
            counts[1][i] = ls2->locus.base_counts_high_qual[i];
        }

        dw->metrics.total++;

        ++dw->pair_stats[pi].total;

        if (! (ls1->is_next && ls2->is_next))
            diff_state = AMBIGUOUS, cacheable = 1, cache_was_set = 1;

        else
            diff_state = 
                cached_dirichlet_diff(counts[0], counts[1], &dw->bep,
                                      &cacheable, &cache_was_set);

        ++dw->pair_stats[pi].dist_count[diff_state];
        dw->pair_stats[pi].cacheable += cacheable;
        dw->pair_stats[pi].cache_was_set += cache_was_set;

        dw->metrics.cacheable += cacheable;
        dw->metrics.cache_was_set += cache_was_set;

        if (diff_state == CHANGED)
        {
            /* Finish sampling and do full distance marginal estimation */
            for (i = 0; i != 2; ++i)
            {
                unsigned perm[] = { 0, 1, 2, 3 };
                struct distrib_points *dst = dw->bep.dist[i];
                /* Generate all points */
                update_points_gen_params(dst, counts[i], perm);
                POINT *p,
                    *pb = dst->points.buf,
                    *pe = dst->points.buf + max_sample_points;
                for (p = pb; p != pe; p += GEN_POINTS_BATCH)
                    dst->pgen.gen_point(dst->pgen.points_gen_par, p);

                dst->points.size = max_sample_points;
                
                /* Generate all weights */
                double *w,
                    *wb = dst->weights.buf,
                    *we = dst->weights.buf + max_sample_points;
                for (w = wb, p = pb; w != we; 
                     w += GEN_POINTS_BATCH, p += GEN_POINTS_BATCH)
                    dst->pgen.weight(p, dst->pgen.points_gen_par, w);

                dst->weights.size = max_sample_points;
            }
            /* Compute weighted square distances (max value of 2).
               e.g.
               minimum simplex distance on [0,1.00] scale is 0.25.
               minimum real    distance on [0,1.41] scale is 0.35355
               minimum sq_real distance on [0,2.00] scale is 0.12499
            */

            compute_wsq_dist(dw->bep.dist[0]->points.buf, 
                             dw->bep.dist[0]->weights.buf,
                             dw->bep.dist[1]->points.buf, 
                             dw->bep.dist[1]->weights.buf,
                             max_sample_points,
                             dw->square_dist_buf,
                             dw->weights_buf);

            double test_quantile = 1.0 - posterior_confidence, test_quantile_value;

            /* Compute the test distance quantile (relative to the
               dist, not squared dist) */
            compute_dist_wquantiles(dw->square_dist_buf,
                                    dw->weights_buf,
                                    max_sample_points,
                                    &test_quantile,
                                    1,
                                    &test_quantile_value);

            if (test_quantile_value > min_dirichlet_dist)
            {
                ++dw->pair_stats[pi].confirmed_changed;
                ls1->confirmed_changed = 1;
                ls2->confirmed_changed = 1;

                if (out_buf)
                {
                    compute_dist_wquantiles(dw->square_dist_buf,
                                            dw->weights_buf,
                                            max_sample_points,
                                            quantiles,
                                            n_quantiles,
                                            dw->dist_quantile_values);            
                    
                    /* quantile values are in squared distance terms, and
                       in the [0, 1] x 4 space of points.  The maximum
                       distance between such points is sqrt(2.0).  We want
                       to re-scale it to be 1. */
                    unsigned q;
                    for (q = 0; q != n_quantiles; ++q)
                        dw->dist_quantile_values[q] *= ONE_OVER_SQRT2;
                    
                    print_distance_quantiles(dw->lslist[gs].locus.reference, 
                                             dw->lslist[gs].locus.position, dw, 
                                             pi, dw->dist_quantile_values, 
                                             out_buf);
                }
            }
        }
        
    }
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

            compute_sq_dist(points1, points2, max_sample_points, n_events,
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
                                               all_events, n_events, dw->lslist, buf);
            }

            delete[] points1;
            delete[] points2;
            delete[] alpha1;
            delete[] alpha2;
        }
        delete[] all_events;
    }
}


// find gs, and init the is_next field for all samples
// this is run
size_t init_global_and_next(locus_sampling *lslist)
{
    size_t gs = 0;
    for (size_t s = 0; s != samples.n; ++s)
    {
        if (lslist[s].current != lslist[s].end
            && (lslist[gs].current == lslist[gs].end
                || cmp_pair_ordering(&lslist[s].locus_ord, 
                                     &lslist[gs].locus_ord) < 0)
            )
            gs = s;
    }
    bool global_at_end = (lslist[gs].current == lslist[gs].end);
    unsigned s;
    for (s = 0; s != samples.n; ++s)
        lslist[s].is_next = 
            lslist[s].current != lslist[s].end
            && (! global_at_end)
            && cmp_pair_ordering(&lslist[s].locus_ord,
                                 &lslist[gs].locus_ord) == 0;

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
void dist_worker(void *par,
                 const struct managed_buf *in_bufs,
                 struct managed_buf *out_bufs)
{
    struct dist_worker_input *dw = (struct dist_worker_input *)par;

    struct timespec worker_start_time;
    clock_gettime(CLOCK_REALTIME, &worker_start_time);
    dw->metrics = { 0, 0, 0 };

    unsigned i = 0;
    struct managed_buf 
        *dist_buf = worker_options.do_dist ? &out_bufs[i++] : NULL,
        *comp_buf = worker_options.do_comp ? &out_bufs[i++] : NULL,
        *indel_buf = worker_options.do_indel ? &out_bufs[i++] : NULL;

    size_t gs = 0, s;
    struct locus_sampling *ls;
    for (s = 0; s != samples.n; ++s)
    {    
        ls = &dw->lslist[s];
        ls->current = in_bufs[s].buf;
        ls->end = in_bufs[s].buf + in_bufs[s].size;
        ((struct points_gen_par *)ls->distp.pgen.points_gen_par)->post_counts = 
            &ls->locus.counts;
    }

    // before main loop, initialize loci and samples
    for (s = 0; s != samples.n; ++s)
    {
        if (dw->lslist[s].current == dw->lslist[s].end) continue;
        update_pileup_locus(&samples.atts[s].nuc_stats, &dw->lslist[s]);
    }

    init_pseudo_locus(&dw->pseudo_sample);

    gs = init_global_and_next(dw->lslist);


    // main loop for computing pairwise distances
    // though slightly inefficient, just construct null_points here
    // !!! here, use worker[0] as a proxy.  this is a design flaw
    // owing to the fact that there are multiple workers used, all with the same
    // basic settings
    locus_sampling null_sd;
    alloc_pileup_locus(&null_sd);

    char null_pileup[] = "chr1\t1\tA\t0\t\t\n";
    null_sd.current = null_pileup;
    null_sd.end = null_pileup + sizeof(null_pileup);
        
    update_pileup_locus(&samples.atts[0].nuc_stats, &null_sd);

    COMP_QV comp_quantile_values;
    double comp_means[NUM_NUCS];

    while (dw->lslist[gs].current != dw->lslist[gs].end)
    {
        /* this will strain the global mutex, but only during the hash
           loading phase. */
        if (! dw->bep.points_hash_frozen)
            dw->bep.points_hash_frozen = points_hash_frozen();

        if (! dw->bep.bounds_hash_frozen)
            dw->bep.bounds_hash_frozen = freeze_bounds_hash();

        if (dist_buf || comp_buf)
            next_distance_quantiles_aux(dw, gs, dist_buf);

        if (indel_buf)
            next_indel_distance_quantiles_aux(dw, gs, indel_buf);

        if (comp_buf)
        {
            for (s = 0; s != samples.n; ++s)
            {
                struct locus_sampling *sam = &dw->lslist[s];
                if (sam->is_next && sam->confirmed_changed)
                {
                    unsigned d;
                    for (d = 0; d != NUM_NUCS; ++d)
                    {
                        compute_marginal_wquantiles((double *)sam->distp.points.buf,
                                                    sam->distp.weights.buf,
                                                    max_sample_points,
                                                    NUM_NUCS,
                                                    d,
                                                    quantiles,
                                                    n_quantiles,
                                                    comp_quantile_values[d]);
                        comp_means[d] = 
                            compute_marginal_mean((double *)sam->distp.points.buf,
                                                  sam->distp.weights.buf,
                                                  max_sample_points,
                                                  NUM_NUCS,
                                                  d);
                    }
                    
                    print_basecomp_quantiles(comp_quantile_values,
                                             comp_means,
                                             n_quantiles,
                                             samples.atts[s].label,
                                             &dw->lslist[s],
                                             comp_buf);
                }
            }
        }
        for (s = 0; s != samples.n; ++s)
            if (dw->lslist[s].is_next) 
                refresh_locus(&samples.atts[s].nuc_stats, &dw->lslist[s]);

        gs = init_global_and_next(dw->lslist);

        /* refresh the pseudo sample */
        update_pseudo_locus(dw->lslist[gs].locus.reference_base, &dw->pseudo_sample);
    }   

    if (dw->do_print_progress)
    {
        struct timespec now;
        clock_gettime(CLOCK_REALTIME, &now);
        unsigned elapsed = now.tv_sec - start_time.tv_sec;

        time_t cal = time(NULL);
        char *ts = strdup(ctime(&cal));
        ts[strlen(ts)-1] = '\0';
        fprintf(stdout, 
                "%s (%02i:%02i:%02i elapsed): Finished processing %s %i\n", 
                ts,
                elapsed / 3600,
                (elapsed % 3600) / 60,
                elapsed % 60,
                dw->lslist[gs].locus.reference,
                dw->lslist[gs].locus.position);
        fflush(stdout);
        free(ts);
    }

    accumulate_pair_stats(dw->pair_stats);

    // struct timespec worker_end_time;
    // clock_gettime(CLOCK_REALTIME, &worker_end_time);
    // unsigned elapsed = worker_end_time.tv_sec - worker_start_time.tv_sec;

    // fprintf(stderr, "%Zu\t%u\t%u\t%u\t%u\n", 
    //         dw->thread_num, elapsed, dw->metrics.total,
    //         dw->metrics.cacheable, dw->metrics.cache_was_set);

    // print_primary_cache_size();
    // print_cache_stats();

    free_pileup_locus(&null_sd);
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
