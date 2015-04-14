#include "tools.h"
#include "usage_strings.h"
#include "defs.h"
#include "dist_worker.h"
#include "pileup_tools.h"
#include "yepLibrary.h"

extern "C" {
#include "dict.h"
#include "range_line_reader.h"
#include "locus.h"
#include "dirichlet_points_gen.h"
}

#define MIN(a,b) ((a) < (b) ? (a) : (b))

int dist_usage()
{
    fprintf(stderr, 
            "\nUsage: dep dist [options] samples.rdb diststats.rdb contig_order.rdb\n" 
            "Options:\n\n"
            "-d STRING   (optional) name of output distance file.  If absent, do not perform distance calculation.\n"
            "-c STRING   (optional) name of output composition file.  If absent, do not output composition marginal estimates.\n"
            "-i STRING   (optional) name of output indel distance file.  If absent, do not perform indel distance calculation.\n"
            "-s STRING   name of sample_pairings.rdb file (see description below).  Required if -d or -i are given.\n"
            "-r STRING   input range file (lines 'contig<tab>start<tab>end') to process [blank] (blank = process whole file)\n"
            "-x STRING   (optional) name of output summary statistics file of differing loci.\n"
            "-q INT      minimum quality score to include bases as evidence [%s]\n"
            "-F STRING   Fastq offset type if known (one of Sanger, Solexa, Illumina13, Illumina15) [None]\n"
            "-g <empty>  if present, print extra pileup fields in the output distance file [absent]\n"
            "-t INT      number of threads to use [1]\n"
            "-m INT      number bytes of memory to use [%Zu]\n"
            // "-v <empty>  if present, be verbose [absent]\n"
            "\n"
            "These options affect how sampling is done.\n"
            "\n"
            "-D STRING   distance quantiles string [\"0.005,0.05,0.5,0.95,0.995\"]\n"
            "-C STRING   composition quantiles string [\"0.005,0.05,0.5,0.95,0.995\"]\n"
            "\n"
            "samples.rdb has lines of <sample_id><tab></path/to/sample.jpd><tab></path/to/sample.pileup>\n"
            "sample_pairings.rdb has lines of <sample_id><tab><sample_id>\n"
            "defining which pairs of samples are to be compared.  <sample_id> correspond\n"
            "with those in samples.rdb\n"
            "\n"
            "diststats.rdb is a file of pre-computed pairwise dirichlet distribution distances, generated\n"
            "by dep diststats\n"
            "\n"
            "contig_order.rdb has lines of <contig><tab><ordering>\n"
            "It must be consistent and complete with the orderings of the contigs mentioned in\n"
            "all pileup input files\n"
            ,
            Usage::quality_string,
            1024 * 1024 * 1024 * 4l // 4 GB memory
            );
    return 1;
}

extern char *optarg;
extern int optind;

int main_dist(int argc, char **argv)
{
    size_t n_threads = 1;
    size_t max_mem = 1024l * 1024l * 1024l * 4l;

    struct posterior_settings pset;
    pset.min_quality_score = 5;

    int print_pileup_fields = 0;

    const char *dist_quantiles_string = "0.005,0.05,0.5,0.95,0.995";
    const char *comp_quantiles_string = "0.005,0.05,0.5,0.95,0.995";

    const char *dist_file = NULL;
    const char *comp_file = NULL;
    const char *indel_dist_file = NULL;
    const char *fastq_type = NULL;
    const char *sample_pairings_file = NULL;
    const char *summary_stats_file = NULL;
    const char *query_range_file = NULL;

    char c;
    while ((c = getopt(argc, argv, "d:c:i:s:x:r:t:m:q:F:D:C:g")) >= 0)
    {
        switch(c)
        {
        case 'd': dist_file = optarg; break;
        case 'c': comp_file = optarg; break;
        case 'i': indel_dist_file = optarg; break;
        case 's': sample_pairings_file = optarg; break;
        case 'x': summary_stats_file = optarg; break;
        case 'r': query_range_file = optarg; break;
        case 't': n_threads = (size_t)atof(optarg); break;
        case 'm': max_mem = (size_t)atof(optarg); break;
        case 'q': pset.min_quality_score = (size_t)atoi(optarg); break;
        case 'F': fastq_type = optarg; break;
        // case 'v': verbose = true; break;
        case 'D': dist_quantiles_string = optarg; break;
        case 'C': comp_quantiles_string = optarg; break;
        case 'g': print_pileup_fields = 1; break;
        default: return dist_usage(); break;
        }
    }
    if (argc - optind != 3) return dist_usage();

    const char *samples_file = argv[optind];
    const char *diststats_file = argv[optind + 1];
    const char *contig_order_file = argv[optind + 2];

    setvbuf(stdout, NULL, _IONBF, 0);
    printf("\n"); /* So progress messages don't interfere with shell prompt. */

    // consistency checks
    if (! sample_pairings_file && (dist_file || indel_dist_file))
    {
        // either provide dist_file and sample_pairings_file, or none at all
        fprintf(stderr, "If -d or -i are given, you must also provide -S.\n");
        exit(5);
    }

    if (! dist_file && ! comp_file && ! indel_dist_file)
    {
        fprintf(stderr, "Error: You must provide at least one of -d or -c or -i.  "
                "Otherwise, there is nothing to calculate\n");
        exit(5);
    }

    FILE *summary_stats_fh = open_if_present(summary_stats_file, "w");
    
    gsl_set_error_handler_off();

    FILE *diststats_fh = fopen(diststats_file, "r");
    if (! diststats_fh)
    {
        fprintf(stderr, "Error: couldn't open diststats file %s for reading.\n", diststats_file);
        exit(1);
    }

    printf("Parsing diststats...");
    unsigned max1_alpha, max2_alpha;
    parse_diststats_header(diststats_fh, &pset, &max1_alpha, &max2_alpha);
    dirichlet_diff_init(max1_alpha, max2_alpha);
    parse_diststats_body(diststats_fh, max1_alpha, max2_alpha);
    fclose(diststats_fh);
    printf("done.\n");

    FILE *samples_fh = open_if_present(samples_file, "r");
    char jpd[1000], pileup[1000], label[1000];
    size_t n_samples = 0, n_alloc = 0;

    struct sample_attributes *sample_attrs_tmp = NULL, *sample_attrs = NULL;
    size_t max_label_len = sizeof((*sample_attrs_tmp).label_string);
    while (! feof(samples_fh))
    {
        (void)fscanf(samples_fh, "%s\t%s\t%s\n", label, jpd, pileup);
        ALLOC_GROW_TYPED(sample_attrs_tmp, n_samples + 1, n_alloc);
        if (strlen(label) > max_label_len)
        {
            fprintf(stderr, "Error: sample label must be less than %Zu characters\n"
                    "Label was \"%s\"\n", max_label_len, label);
            exit(1);
        }
        init_sample_attributes(jpd, label, pileup, &sample_attrs_tmp[n_samples]);
        ++n_samples;
    }
    fclose(samples_fh);

    int offset;
    if (fastq_type) offset = fastq_type_to_offset(fastq_type);
    else
    {
        size_t sz = 1024 * 1024 * 64;
        char *buf = (char *)malloc(sz);

        /* using 'pileup' here.  will be the name of the last pileup
           file mentioned in the samples file.  */
        offset = fastq_offset(pileup, buf, sz);
        free(buf);
    }

    if (offset == -1)
    {
        fprintf(stderr, "Could not determine fastq type of this pileup file.\n");
        return 1;
    }

    PileupSummary::set_offset(offset);

    FILE *dist_fh = open_if_present(dist_file, "w");
    FILE *comp_fh = open_if_present(comp_file, "w");
    FILE *indel_fh = open_if_present(indel_dist_file, "w");


    /* hack to initialize both dist and comp quantiles */
    struct {
        const char *quant_string;
        double *quant;
        size_t *n_quant;
    } Q[] = { 
        { dist_quantiles_string, pset.dist_quantiles, &pset.n_dist_quantiles }, 
        { comp_quantiles_string, pset.comp_quantiles, &pset.n_comp_quantiles }
    };

    size_t rval = 0;
    double qval;
    int i, pos, off;
    pset.n_dist_quantiles = 0;
    pset.n_comp_quantiles = 0;
    for (i = 0; i != 2; ++i)
    {
        pos = 0;
        off = 0;
        while (Q[i].quant_string[off] != '\0' &&
               (rval = sscanf(Q[i].quant_string + off, "%lf%n", &qval, &pos)) == 1)
        {
            if (*Q[i].n_quant == MAX_NUM_QUANTILES) break;
            Q[i].quant[(*Q[i].n_quant)++] = qval;
            off += pos + (Q[i].quant_string[off + pos] == ',' ? 1 : 0);
        }
        if (rval != 1 || *Q[i].n_quant == 0 
            || *Q[i].n_quant == MAX_NUM_QUANTILES)
        {
            fprintf(stderr, "Error: something amiss with quantiles\n%s\n",
                    Q[i].quant_string);
            exit(1);
        }
    }

    // parse pairings file
    size_t n_pairings, *pair_sample1, *pair_sample2;
    {
        int *sample_remap = (int *)malloc(n_samples * sizeof(int));
        size_t p, s, n = 0, na = 0;
        for (s = 0; s != n_samples; ++s) sample_remap[s] = -1;

        FILE *sample_pairings_fh = open_if_present(sample_pairings_file, "r");
        std::vector<size_t> pair1, pair2;
        std::vector<size_t> *trg[] = { & pair1, & pair2 };
        char label[2][100];
        while (! feof(sample_pairings_fh))
        {
            (void)fscanf(sample_pairings_fh, "%s\t%s\n", label[0], label[1]);
            for (p = 0; p != 2; ++p)
            {
                for (s = 0; s != n_samples; ++s)
                    if (strcmp(sample_attrs_tmp[s].label_string, label[p]) == 0)
                    {
                        if (sample_remap[s] == -1)
                        {
                            sample_remap[s] = n;
                            ALLOC_GROW_TYPED(sample_attrs, n + 1, na);
                            sample_attrs[n] = sample_attrs_tmp[s];
                            ++n;
                        }
                        (*trg[p]).push_back(sample_remap[s]);
                        break;
                    }
                if (s == n_samples)
                {
                    fprintf(stderr, "Error: label %s in %s not found "
                            "in samples file %s\n",
                            label[p], sample_pairings_file, samples_file);
                    exit(10);
                }
            }
        }
        fclose(sample_pairings_fh);
        free(sample_remap);
        
        n_pairings = pair1.size();
        pair_sample1 = (size_t *)malloc(n_pairings * sizeof(size_t));
        pair_sample2 = (size_t *)malloc(n_pairings * sizeof(size_t));
        std::copy(pair1.begin(), pair1.end(), pair_sample1);
        std::copy(pair2.begin(), pair2.end(), pair_sample2);
        n_samples = n;
    }
    free(sample_attrs_tmp);

    /* */
    
    /* 0. parse contig order file */
    char contig[1024];
    unsigned index;
    FILE *contig_order_fh = open_if_present(contig_order_file, "r");
    while (! feof(contig_order_fh))
    {
        int n = fscanf(contig_order_fh, "%s\t%u\n", contig, &index);
        if (n != 2)
        {
            fprintf(stderr, "Error: contig order file %s doesn't have the proper format\n",
                    contig_order_file);
            exit(1);
        }
        dict_add_item(contig, index);
    }
    fclose(contig_order_fh);
    
    dict_build();

    size_t t, s; /* threads, samples */

    struct dist_worker_input *worker_buf = (struct dist_worker_input *)
        malloc(n_threads * sizeof(struct dist_worker_input));
    
    void **worker_inputs = (void **)malloc(n_threads * sizeof(void *));
    // initialize all generic structures in the worker_inputs

    /* This adjustment makes max_sample_points a multiple of GEN_POINTS_BATCH */
    pset.max_sample_points += GEN_POINTS_BATCH - (pset.max_sample_points % GEN_POINTS_BATCH);

    for (t = 0; t != n_threads; ++t)
    {
        worker_buf[t].sample_atts = sample_attrs;
        worker_buf[t].bep.max1 = max1_alpha;
        worker_buf[t].bep.max2 = max2_alpha;
        worker_buf[t].bep.batch_size = GEN_POINTS_BATCH;
        worker_buf[t].bep.pset = &pset; 
        worker_buf[t].thread_num = t;
        worker_buf[t].n_samples = n_samples;
        worker_buf[t].n_sample_pairs = n_pairings;
        worker_buf[t].randgen = gsl_rng_alloc(gsl_rng_taus);
        worker_buf[t].square_dist_buf = 
            (double *)malloc(sizeof(double) * pset.max_sample_points);
        worker_buf[t].weights_buf =
            (double *)malloc(sizeof(double) * pset.max_sample_points);
        worker_buf[t].pair_stats = 
            (struct pair_dist_stats *)calloc(n_pairings, sizeof(struct pair_dist_stats));
        worker_buf[t].print_pileup_fields = print_pileup_fields;
        worker_buf[t].do_dist = (dist_fh != NULL);
        worker_buf[t].do_comp = (comp_fh != NULL);
        worker_buf[t].do_indel = (indel_fh != NULL);
        worker_buf[t].pair_sample1 = pair_sample1;
        worker_buf[t].pair_sample2 = pair_sample2;
    }

    for (t = 0; t != n_threads; ++t)
        worker_inputs[t] = &worker_buf[t];

    printf("Creating in-memory index of input files...");
    size_t scan_thresh_size = 1e6;
    file_bsearch_init(init_locus, scan_thresh_size);
    struct file_bsearch_index *ix = (struct file_bsearch_index *)
        malloc(n_samples * sizeof(struct file_bsearch_index));

    for (s = 0; s != n_samples; ++s)
        ix[s] = file_bsearch_make_index(sample_attrs[s].fh);
    printf("done.\n");

    /* initialize queries */
    struct pair_ordering_range *queries, *q, *qend;
    size_t n_queries;

    if (query_range_file)
        queries = parse_query_ranges(query_range_file, &n_queries);
    else
    {
        /* simply set the 'query' to the entire file */
        n_queries = 1;
        queries = (struct pair_ordering_range *)
            malloc(sizeof(struct pair_ordering_range));
        queries[0] = { ix[0].root->span_beg, ix[0].root->span_end };
        queries[0].end.lo--; /* necessary for this pseudo-query to fit in the index root */
    }
    q = queries;
    qend = queries + n_queries;

    struct range_line_reader_par reader_par = {
        ix, n_samples, q, qend, init_locus, 1
    };

    struct dist_worker_offload_par offload_par = {
        dist_fh, comp_fh, indel_fh
    };

    size_t n_output_files = 
        (dist_fh ? 1 : 0)
        + (comp_fh ? 1 : 0)
        + (indel_fh ? 1 : 0);

    /* initialize beta quantile estimation procedure */
    printf("Precomputing confidence interval statistics...");
    init_beta(pset.beta_confidence, n_threads);
    printf("done.\n");

    /* this should be revisited if it turns out that threads are
       stalling */
    size_t n_extra = n_threads;
    struct thread_queue *tqueue =
        thread_queue_init(range_line_reader, &reader_par,
                          dist_worker, worker_inputs,
                          dist_offload, &offload_par,
                          n_threads,
                          n_extra,
                          n_samples,
                          n_output_files,
                          max_mem);

    enum YepStatus status = yepLibrary_Init();
    assert(status == YepStatusOk);

    printf("Starting input processing.\n");
    thread_queue_run(tqueue);
    thread_queue_free(tqueue);

    /* Aggregate total pair statistics */
    struct pair_dist_stats *all_pair_stats = 
        (struct pair_dist_stats *)calloc(n_pairings, sizeof(struct pair_dist_stats));

    unsigned n_states = sizeof(all_pair_stats[0].dist_count) / sizeof(all_pair_stats[0].dist_count[0]);

    if (summary_stats_fh)
    {
        unsigned p;
        fprintf(summary_stats_fh, "sample1\tsample2");
        fprintf(summary_stats_fh, "\t%s\t%s\t%s", "cache_hit", "cache_miss", "confirmed_changed");
        for (s = 0; s != n_states; ++s)
            fprintf(summary_stats_fh, "\t%s", fuzzy_state_strings[s]);
        fprintf(summary_stats_fh, "\n");

        for (p = 0; p != n_pairings; ++p)
        {
            for (t = 0; t != n_threads; ++t)
            {
                all_pair_stats[p].confirmed_changed 
                    += worker_buf[t].pair_stats[p].confirmed_changed;
                all_pair_stats[p].cache_hit
                    += worker_buf[t].pair_stats[p].cache_hit;
                all_pair_stats[p].cache_miss
                    += worker_buf[t].pair_stats[p].cache_miss;

                for (s = 0; s != n_states; ++s)
                    all_pair_stats[p].dist_count[s] 
                        += worker_buf[t].pair_stats[p].dist_count[s];
            }

            /* Print out statistics */
            fprintf(summary_stats_fh, "%s\t%s", 
                    sample_attrs[pair_sample1[p]].label_string,
                    sample_attrs[pair_sample2[p]].label_string);
            fprintf(summary_stats_fh, "\t%zu", all_pair_stats[p].cache_hit);
            fprintf(summary_stats_fh, "\t%zu", all_pair_stats[p].cache_miss);
            fprintf(summary_stats_fh, "\t%zu", all_pair_stats[p].confirmed_changed);
            for (s = 0; s != n_states; ++s)
                fprintf(summary_stats_fh, "\t%zu", all_pair_stats[p].dist_count[s]);
            fprintf(summary_stats_fh, "\n");

        }
        fclose(summary_stats_fh);
    }

    free(ix);
    for (t = 0; t != n_threads; ++t)
    {
        gsl_rng_free(worker_buf[t].randgen);
        free(worker_buf[t].pair_stats);
        free(worker_buf[t].square_dist_buf);
        free(worker_buf[t].weights_buf);
    }

    free(all_pair_stats);
    free(worker_buf);
    free(worker_inputs);

    if (dist_fh) fclose(dist_fh);
    if (comp_fh) fclose(comp_fh);
    if (indel_fh) fclose(indel_fh);

    for (s = 0; s != n_samples; ++s) fclose(sample_attrs[s].fh);

    free(pair_sample1);
    free(pair_sample2);
    free(sample_attrs);
    dirichlet_diff_free();
    printf("Finished.\n");

    return 0;
}
