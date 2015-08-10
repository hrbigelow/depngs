#include <string.h>

#include "common_tools.h"
#include "defs.h"
#include "locus_diff.h"
#include "bam_sample_info.h"

#include "locus.h"
#include "dirichlet_points_gen.h"


#define MIN(a,b) ((a) < (b) ? (a) : (b))

int dist_usage()
{
    fprintf(stderr, 
            "\n\n"
            "Usage: dep dist [options] samples.rdb sample_pairings.rdb ref.fasta\n" 
            "Options:\n\n"
            "FILES\n"
            "-d STRING   name of output distance file.  If absent, skip distance calculation.\n"
            "-c STRING   name of output composition file.  If absent, skip composition estimates.\n"
            "-i STRING   name of output indel distance file.  If absent, skip indel distance calculation.\n"
            "-l STRING   input locus range file (lines 'contig<tab>start<tab>end') to process.\n"
            "              If absent, process all input.\n"
            "-x STRING   name of output summary statistics file of differing loci.\n"
            "\n"
            "STATISTICAL PARAMETERS\n"
            "-y FLOAT    MIN_DIST, min mutation distance (0-1 scale) to call as changed [0.2]\n"
            "-X FLOAT    POST_CONF, confidence for -y [0.99]\n"
            "-Z FLOAT    BETA_CONF, confidence for binomial estimation [0.9999]\n"
            "-f INT      MAX_POINTS, max # of sample point pairs for binomial test [10000]\n"
            "-p FLOAT    PRIOR_ALPHA, alpha value for dirichlet prior [0.1]\n"
            "-Q STRING   dist/comp/indel quantiles to report [\"0.005,0.05,0.5,0.95,0.995\"]\n"
            "\n"
            "OTHER\n"
            "-q INT      minimum quality score to include bases as evidence [5]\n"
            "-g <empty>  if present, print extra pileup fields in the output distance file [absent]\n"
            "-t INT      number of threads to use [1]\n"
            "-R INT      number of readers to use [same as -t] (See NOTE)\n"
            "-m INT      number bytes of memory to use [10e9]\n"
            "\n"
            "\n"
            "samples.rdb has lines of <sample_id><tab></path/to/sample.bam>\n"
            "\n"
            "sample_pairings.rdb has lines of <sample_id><tab><sample_id>\n"
            "defining which pairs of samples are to be compared.  <sample_id>\n"
            "correspond with those in samples.rdb.  The special\n"
            "<sample_id> \"REF\" may be supplied as the second sample in the\n"
            "pair.  This indicates to do comparisons with a conceptual sample\n"
            "that is identical to the reference base at every locus.\n"
            "\n"
            "ref.fasta is the fasta-formatted reference genome.  all bam file\n"
            "inputs must be aligned to this reference.  there must be a fasta\n"
            "index file (produced by samtools faidx) called <ref.fasta>.fai\n"
            "\n"
            "NOTE: -R may be used to restrict how many threads are allowed to\n"
            "read input concurrently, which may improve performance if there\n"
            "are many threads and reading is a bottleneck\n"
            "\n"
            );
    return 1;
}

extern char *optarg;
extern int optind;


int
main_dist(int argc, char **argv)
{
    size_t n_threads = 1;
    size_t max_mem = 10e9;
    unsigned n_readers = n_threads;
    int n_readers_set = 0;

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
    const char *summary_stats_file = NULL;
    const char *query_range_file = NULL;

    double prior_alpha = 0.1;
#define SQRT2 1.41421356237309504880

    char c;
    while ((c = getopt(argc, argv, "d:c:i:l:x:t:R:f:y:X:Z:p:t:m:q:Q:g")) >= 0) {
        switch(c) {
        case 'd': dist_file = optarg; break;
        case 'c': comp_file = optarg; break;
        case 'i': indel_dist_file = optarg; break;
        case 'l': query_range_file = optarg; break;
        case 'x': summary_stats_file = optarg; break;

        case 'f': max_sample_points = 
                (unsigned)strtod_errmsg(optarg, "-f (max_sample_points)"); break;
        case 'y': min_dirichlet_dist = 
                SQRT2 * strtod_errmsg(optarg, "-y (min_dirichlet_dist)"); break;
        case 'X': post_confidence = strtod_errmsg(optarg, "-X (post_confidence)"); break;
        case 'Z': beta_confidence = strtod_errmsg(optarg, "-Z (beta_confidence)"); break;
        case 'p': prior_alpha = strtod_errmsg(optarg, "-p (prior_alpha)"); break;

        case 't': n_threads = strtol_errmsg(optarg, "-t (n_threads)"); break;
        case 'R': 
            n_readers = strtol_errmsg(optarg, "-R (n_readers)"); 
            n_readers_set = 1;
            break;
        case 'm': max_mem = (size_t)strtod_errmsg(optarg, "-m (max_mem)"); break;
        case 'q': min_quality_score = strtol_errmsg(optarg, "-q (min_quality_score)"); break;
        case 'Q': quantiles_string = optarg; break;
        case 'g': print_pileup_fields = 1; break;
        default: return dist_usage(); break;
        }
    }
    if (argc - optind != 3) return dist_usage();

    if (! n_readers_set)
        n_readers = n_threads;

    /* This adjustment makes max_sample_points a multiple of GEN_POINTS_BATCH */
    max_sample_points += GEN_POINTS_BATCH - (max_sample_points % GEN_POINTS_BATCH);

    const char *samples_file = argv[optind];
    const char *sample_pairs_file = argv[optind + 1];
    const char *fasta_file = argv[optind + 2];

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

#define BYTES_PER_POINT sizeof(double) * NUM_NUCS

    /* this is just an empirically based estimate */
#define FRAGMENTATION_FACTOR 0.85

    /* */
#define INPUT_MEM_FRACTION 0.4

    /* allot a fraction of total memory to input buffers. */
    unsigned long max_input_mem = max_mem * INPUT_MEM_FRACTION;
    unsigned long max_dir_cache_items = 
        FRAGMENTATION_FACTOR * (max_mem - max_input_mem) 
        / (max_sample_points * BYTES_PER_POINT);
    unsigned long max_bounds_cache_items = 5000;

    locus_diff_init(post_confidence, beta_confidence, 
                    min_dirichlet_dist, max_sample_points,
                    max_dir_cache_items, max_bounds_cache_items,
                    n_threads, prior_alpha,
                    samples_file, sample_pairs_file, fasta_file,
                    min_quality_score, quantiles_string,
                    (dist_fh != NULL), (comp_fh != NULL), (indel_fh != NULL),
                    print_pileup_fields);


    /* initialize beta quantile estimation procedure */

    set_points_hash_flag(1);

    printf("Precomputing difference hash...");
    prepopulate_bounds_keys(n_threads);
    printf("done.\n");

    set_points_hash_flag(0);

    struct thread_queue *tqueue =
        locus_diff_tq_init(query_range_file, 
                           fasta_file,
                           n_threads, n_readers, max_input_mem,
                           dist_fh, comp_fh, indel_fh);

    printf("Starting input processing.\n");
    thread_queue_run(tqueue);
    thread_queue_free(tqueue);

    if (summary_stats_file)
        print_pair_stats(summary_stats_file);

    locus_diff_tq_free();

    if (dist_fh) fclose(dist_fh);
    if (comp_fh) fclose(comp_fh);
    if (indel_fh) fclose(indel_fh);

    locus_diff_free();

    printf("Finished.\n");

    return 0;
}
