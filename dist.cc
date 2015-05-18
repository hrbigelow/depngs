#include "tools.h"
#include "usage_strings.h"
#include "defs.h"
#include "dist_worker.h"
#include "pileup_tools.h"

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
            "\nUsage: dep dist [options] samples.rdb sample_pairings.rdb contig_order.rdb\n" 
            "Options:\n\n"
            "FILES\n"
            "-d STRING   (optional) name of output distance file.  If absent, do not perform distance calculation.\n"
            "-c STRING   (optional) name of output composition file.  If absent, do not output composition marginal estimates.\n"
            "-i STRING   (optional) name of output indel distance file.  If absent, do not perform indel distance calculation.\n"
            "-r STRING   input range file (lines 'contig<tab>start<tab>end') to process [blank] (blank = process whole file)\n"
            "-x STRING   (optional) name of output summary statistics file of differing loci.\n"
            "\n"
            "STATISTICAL PARAMETERS\n"
            "-y FLOAT    MIN_DIST, min mutation distance (0-1 scale) to call as changed [0.2]\n"
            "-X FLOAT    POST_CONF, confidence for -y [0.99]\n"
            "-Z FLOAT    BETA_CONF, confidence for binomial estimation [0.9999]\n"
            "-f INT      MAX_POINTS, max # of sample point pairs for binomial test [10000]\n"
            "-p FLOAT    PRIOR_ALPHA, alpha value for dirichlet prior [0.1]\n"
            "\n"
            "OTHER\n"
            "-q INT      minimum quality score to include bases as evidence [5]\n"
            "-F STRING   Fastq offset type if known (one of Sanger, Solexa, Illumina13, Illumina15) [None]\n"
            "-g <empty>  if present, print extra pileup fields in the output distance file [absent]\n"
            "-t INT      number of threads to use [1]\n"
            "-R INT      number of readers to use [2] (See NOTE)\n"
            "-m INT      number bytes of memory to use [50e9]\n"
            // "-v <empty>  if present, be verbose [absent]\n"
            "\n"
            "These options affect how sampling is done.\n"
            "\n"
            "-Q STRING   quantiles string [\"0.005,0.05,0.5,0.95,0.995\"]\n"
            "\n"
            "samples.rdb has lines of <sample_id><tab></path/to/sample.jpd><tab></path/to/sample.pileup>\n"
            "sample_pairings.rdb has lines of <sample_id><tab><sample_id>\n"
            "defining which pairs of samples are to be compared.  <sample_id> correspond\n"
            "with those in samples.rdb.  The special <sample_id> \"PSEUDO\" may be"
            "supplied as the second sample in the pair.  This indicates to do comparisons\n"
            "with a conceptual sample that is identical to the reference base at every locus.\n"
            "\n"
            "contig_order.rdb has lines of <contig><tab><ordering>\n"
            "It must be consistent and complete with the orderings of the contigs mentioned in\n"
            "all pileup input files\n"
            "\n"
            "On machines where 2 or more concurrent reads achieve higher\n"
            "throughput than one read, set -R (number of readers) accordingly\n"
            );
    return 1;
}

extern char *optarg;
extern int optind;

int main_dist(int argc, char **argv)
{
    size_t n_threads = 1;
    size_t max_mem = 50e9;
    unsigned n_readers = 2;

    unsigned min_quality_score = 5;
    double post_confidence = 0.99;
    double beta_confidence = 0.9999;
    double min_dirichlet_dist = 0.2;
    unsigned max_sample_points = 10000;
    int print_pileup_fields = 0;

    const char *quantiles_string = "0.005,0.05,0.5,0.95,0.995";

    const char *dist_file = NULL;
    const char *comp_file = NULL;
    const char *indel_dist_file = NULL;
    const char *fastq_type = NULL;
    const char *summary_stats_file = NULL;
    const char *query_range_file = NULL;

    double prior_alpha = 0.1;

    char c;
    while ((c = getopt(argc, argv, "d:c:i:r:x:t:R:f:y:X:Z:p:t:m:q:F:Q:g")) >= 0)
    {
        switch(c)
        {
        case 'd': dist_file = optarg; break;
        case 'c': comp_file = optarg; break;
        case 'i': indel_dist_file = optarg; break;
        case 'r': query_range_file = optarg; break;
        case 'x': summary_stats_file = optarg; break;

        case 'f': max_sample_points = (size_t)atof(optarg); break;
        case 'y': min_dirichlet_dist = sqrt(2.0) * atof(optarg); break;
        case 'X': post_confidence = atof(optarg); break;
        case 'Z': beta_confidence = atof(optarg); break;
        case 'p': prior_alpha = atof(optarg); break;

        case 't': n_threads = (size_t)atof(optarg); break;
        case 'R': n_readers = (unsigned)atof(optarg); break;
        case 'm': max_mem = (size_t)atof(optarg); break;
        case 'q': min_quality_score = (size_t)atoi(optarg); break;
        case 'F': fastq_type = optarg; break;
        case 'Q': quantiles_string = optarg; break;
        case 'g': print_pileup_fields = 1; break;
        default: return dist_usage(); break;
        }
    }
    if (argc - optind != 3) return dist_usage();

    pileup_init(min_quality_score);

    /* This adjustment makes max_sample_points a multiple of GEN_POINTS_BATCH */
    max_sample_points += GEN_POINTS_BATCH - (max_sample_points % GEN_POINTS_BATCH);

    const char *samples_file = argv[optind];
    const char *sample_pairs_file = argv[optind + 1];
    const char *contig_order_file = argv[optind + 2];

    setvbuf(stdout, NULL, _IONBF, 0);
    printf("\n"); /* So progress messages don't interfere with shell prompt. */

    if (! dist_file && ! comp_file && ! indel_dist_file)
    {
        fprintf(stderr, "Error: You must provide at least one of -d or -c or -i.  "
                "Otherwise, there is nothing to calculate\n");
        exit(5);
    }

    gsl_set_error_handler_off();

    FILE *dist_fh = open_if_present(dist_file, "w");
    FILE *comp_fh = open_if_present(comp_file, "w");
    FILE *indel_fh = open_if_present(indel_dist_file, "w");

    dist_worker_init(post_confidence, min_dirichlet_dist, max_sample_points,
                     samples_file, sample_pairs_file, fastq_type,
                     quantiles_string,
                     (dist_fh != NULL), (comp_fh != NULL), (indel_fh != NULL),
                     print_pileup_fields);

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

#define BYTES_PER_POINT sizeof(double) * NUM_NUCS

    /* this is just an empirically based estimate */
#define FRAGMENTATION_FACTOR 0.85

    /* */
#define INPUT_MEM_FRACTION 0.4

    /* allot 10% of total memory to input buffers. */
    unsigned long max_input_mem = max_mem * INPUT_MEM_FRACTION;
    unsigned long max_dir_cache_items = 
        FRAGMENTATION_FACTOR * (max_mem - max_input_mem) 
        / (max_sample_points * BYTES_PER_POINT);
    unsigned long max_bounds_cache_items = 5000;

    init_dirichlet_points_gen(prior_alpha);

    dirichlet_diff_init(PSEUDO_DEPTH,
                        GEN_POINTS_BATCH,
                        post_confidence, beta_confidence,
                        min_dirichlet_dist,
                        max_sample_points,
                        max_dir_cache_items, 
                        max_bounds_cache_items,
                        n_threads);

    /* initialize beta quantile estimation procedure */
    printf("Precomputing confidence interval statistics...");
    binomial_est_init(beta_confidence, GEN_POINTS_BATCH, 
                      max_sample_points, n_threads);

    printf("done.\n");

    printf("Prepopulating Difference hash...");
    prepopulate_bounds_keys(n_threads);

    struct thread_queue *tqueue =
        dist_worker_tq_init(query_range_file, 
                            n_threads, n_readers, max_input_mem,
                            dist_fh, comp_fh, indel_fh);

    printf("Starting input processing.\n");
    thread_queue_run(tqueue);
    thread_queue_free(tqueue);

    print_pair_stats(summary_stats_file);

    dist_worker_tq_free();

    if (dist_fh) fclose(dist_fh);
    if (comp_fh) fclose(comp_fh);
    if (indel_fh) fclose(indel_fh);

    dist_worker_free();
    dirichlet_diff_free();
    binomial_est_free();

    printf("Finished.\n");

    return 0;
}
