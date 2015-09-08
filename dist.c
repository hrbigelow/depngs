#include <string.h>

#include "common_tools.h"
#include "defs.h"
#include "locus_diff.h"
#include "bam_sample_info.h"
#include "dir_cache.h"
#include "dir_points_gen.h"
#include "dir_diff_cache.h"
#include "timer.h"

static struct {
    struct bam_filter_params bf_par;
    struct binomial_est_params be_par;
    struct dir_cache_params dc_par;
    struct dirichlet_diff_params dd_par;
    struct locus_diff_params ld_par;
    unsigned long max_mem;
    unsigned n_threads;
    unsigned n_max_reading;
} opts = { 
    .bf_par = { 
        .min_base_quality = 5,
        .min_map_quality = 10,
        .rflag_require = 0,
        .rflag_filter = 0
    }, 
    .be_par = {
        .max_sample_points = 10000,
        .post_confidence = 0.99,
        .beta_confidence = 0.9999,
        .min_dirichlet_dist = 0.2,
        .batch_size = GEN_POINTS_BATCH
    },
    .dc_par = {
        .n_bounds = 1e8,
        .min_ct_keep_bound = 3,
        .fasta_file = NULL,
        .n_max_survey_loci = 1e8
    },
    .dd_par = {
        .pseudo_depth = 1e6,
        .prior_alpha = 0.1,
        .xmax = 10000,
        .mode_batch_size = 512,
        .max_bernoulli_trials = 50000
    },
    .ld_par = {
        .do_print_pileup = 0,
        .indel_prior_alpha = 0.1
    },
    .max_mem = 10e9, 
    .n_threads = 1,
    .n_max_reading = 1
};


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
            "\n");

    fprintf(stderr,
            "STATISTICAL PARAMETERS\n"
            "-y FLOAT    MIN_DIST, min mutation distance (0-1 scale) to call as changed [%0.3f]\n"
            "-X FLOAT    POST_CONF, confidence for -y [%0.3f]\n"
            "-Z FLOAT    BETA_CONF, confidence for binomial estimation [%0.5g]\n"
            "-P INT      MAX_POINTS, max # of sample point pairs for binomial test [%d]\n"
            "-p FLOAT    PRIOR_ALPHA, alpha value for dirichlet prior [%0.5g]\n"
            "-C STRING   dist/comp/indel quantiles to report [\"0.005,0.05,0.5,0.95,0.995\"]\n"
            "-S INT      MAX_SURVEY_LOCI, number of loci to survey to build dirichlet statistics[%ld]\n"
            "\n",
            opts.be_par.min_dirichlet_dist,
            opts.be_par.post_confidence,
            opts.be_par.beta_confidence,
            opts.be_par.max_sample_points,
            opts.ld_par.prior_alpha,
            opts.dc_par.n_max_survey_loci
            );

    fprintf(stderr,
            "READ LEVEL FILTERING\n"
            "-Q INT      minimum base quality score to include bases as evidence [%d]\n"
            "-q INT      minimum mapping quality to include a BAM record [%d]\n"
            "-f INT      bits required to be set in BAM flag for inclusion [%d]\n"
            "-F INT      bits required to be unset in BAM flag for inclusion [%d]\n"
            "-R STRING   file listing read groups to require. (absent: do not filter by readgroup) [absent]\n"
            "\n",
            opts.bf_par.min_base_quality,
            opts.bf_par.min_map_quality,
            opts.bf_par.rflag_require,
            opts.bf_par.rflag_filter
            );

    fprintf(stderr,
            "GENERAL\n"
            "-t INT      number of threads to use [%u]\n"
            "-r INT      max number of threads allowed to read input concurrently [same as -t](See NOTE)\n"
            "-m INT      number bytes of memory to use [%lu]\n"
            "-g FLAG     if present, print extra pileup fields in the output distance file [absent]\n"
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
            "NOTE: -r may be used to restrict how many threads are allowed to\n"
            "read input concurrently, which may improve performance if there\n"
            "are many threads and reading is a bottleneck\n"
            "\n",
            opts.n_threads,
            opts.max_mem
);

    return 1;
}

extern char *optarg;
extern int optind;


int
main_dist(int argc, char **argv)
{
    timer_init();
    
    int n_max_reading_set = 0;

    const char *quantiles_string = "0.005,0.05,0.5,0.95,0.995";

    const char *dist_file = NULL;
    const char *comp_file = NULL;
    const char *indel_dist_file = NULL;
    const char *summary_stats_file = NULL;
    const char *query_range_file = NULL;
    const char *readgroup_file = NULL;

#define SQRT2 1.41421356237309504880

    char c;
    while ((c = getopt(argc, argv, "d:c:i:l:x:y:X:Z:P:S:p:C:Q:q:f:F:R:t:r:m:g")) >= 0) {
        switch(c) {
            /* files */
        case 'd': dist_file = optarg; break;
        case 'c': comp_file = optarg; break;
        case 'i': indel_dist_file = optarg; break;
        case 'l': query_range_file = optarg; break;
        case 'x': summary_stats_file = optarg; break;

            /* statistical parameters */
        case 'y': opts.be_par.min_dirichlet_dist = 
                SQRT2 * strtod_errmsg(optarg, "-y (min_dirichlet_dist)"); break;
        case 'X': opts.be_par.post_confidence = strtod_errmsg(optarg, "-X (post_confidence)"); break;
        case 'Z': opts.be_par.beta_confidence = strtod_errmsg(optarg, "-Z (beta_confidence)"); break;
        case 'P': opts.be_par.max_sample_points = 
                (unsigned)strtod_errmsg(optarg, "-f (max_sample_points)"); 
            break;
        case 'p': opts.ld_par.prior_alpha = strtod_errmsg(optarg, "-p (prior_alpha)"); break;
        case 'C': quantiles_string = optarg; break;
        case 'S': opts.dc_par.n_max_survey_loci = strtod_errmsg(optarg, "-S (n_max_survey_loci)"); break;

            /* read-level filtering */
        case 'Q': opts.bf_par.min_base_quality = strtol_errmsg(optarg, "-Q (min_base_quality)"); break;
        case 'q': opts.bf_par.min_map_quality = strtol_errmsg(optarg, "-q (min_map_quality)"); break;
        case 'f': opts.bf_par.rflag_require = strtol_errmsg(optarg, "-f (rflag_require)"); break;
        case 'F': opts.bf_par.rflag_filter = strtol_errmsg(optarg, "-F (rflag_filter)"); break;
        case 'R': readgroup_file = optarg; break;

            /* general */
        case 't': opts.n_threads = strtol_errmsg(optarg, "-t (n_threads)"); break;
        case 'r': 
            opts.n_max_reading = strtol_errmsg(optarg, "-r (n_max_reading)"); 
            n_max_reading_set = 1;
            break;
        case 'm': opts.max_mem = (size_t)strtod_errmsg(optarg, "-m (max_mem)"); break;
        case 'g': opts.ld_par.do_print_pileup = 1; break;

        default: return dist_usage(); break;
        }
    }
    if (argc - optind != 3) return dist_usage();

    if (! n_max_reading_set)
        opts.n_max_reading = opts.n_threads;

    /* This adjustment makes max_sample_points a multiple of GEN_POINTS_BATCH */
    opts.be_par.max_sample_points += 
        GEN_POINTS_BATCH - (opts.be_par.max_sample_points % GEN_POINTS_BATCH);

    const char *samples_file = argv[optind];
    const char *sample_pairs_file = argv[optind + 1];
    const char *fasta_file = argv[optind + 2];

    setvbuf(stdout, NULL, _IONBF, 0);
    printf("\n"); /* So progress messages don't interfere with shell prompt. */

    if (! dist_file && ! comp_file && ! indel_dist_file) {
        fprintf(stderr, "Error: You must provide at least one of -d or -c or -i.  "
                "Otherwise, there is nothing to calculate\n");
        exit(5);
    }

    /* parse readgroups file */
    FILE *readgroup_fh = open_if_present(readgroup_file, "r");
    if (readgroup_fh) {
        opts.bf_par.readgroup_include_hash = 
            init_readgroup_file(readgroup_fh);
        fclose(readgroup_fh);
    }

    gsl_set_error_handler_off();

    FILE *dist_fh = open_if_present(dist_file, "w");
    FILE *comp_fh = open_if_present(comp_file, "w");
    FILE *indel_fh = open_if_present(indel_dist_file, "w");

    opts.ld_par.do_dist = (dist_fh != NULL);
    opts.ld_par.do_comp = (comp_fh != NULL);
    opts.ld_par.do_indel = (indel_fh != NULL);

    /* resolve overlap between parameter sets */
    opts.ld_par.max_sample_points = opts.be_par.max_sample_points;
    opts.ld_par.post_confidence = opts.be_par.post_confidence;
    opts.ld_par.min_dirichlet_dist = opts.be_par.min_dirichlet_dist;
    opts.ld_par.prior_alpha = opts.dd_par.prior_alpha;
    opts.dc_par.fasta_file = fasta_file;
    opts.dc_par.max_sample_points = opts.be_par.max_sample_points;

    /* allot fractions of main memory to points and input buffers.
       bounds and output buffers will be negligible */
#define FRAC_MEM_POINTSETS 0.9
    size_t max_point_sets = opts.max_mem 
        / (sizeof(POINT) * opts.be_par.max_sample_points);
    opts.dc_par.n_point_sets = max_point_sets * FRAC_MEM_POINTSETS;

#define INPUT_MEM_FRACTION 0.05

    /* allot a fraction of total memory to input buffers. */
    unsigned long max_input_mem = opts.max_mem * INPUT_MEM_FRACTION;

    parse_csv_line(quantiles_string, 
                   opts.ld_par.quantiles, 
                   &opts.ld_par.n_quantiles, 
                   MAX_NUM_QUANTILES);

    struct thread_queue *tqueue =
        locus_diff_init(samples_file, sample_pairs_file, 
                        query_range_file, fasta_file,
                        opts.n_threads, opts.n_max_reading, max_input_mem,
                        opts.ld_par, opts.dd_par, opts.be_par, opts.dc_par, opts.bf_par,
                        dist_fh, comp_fh, indel_fh);

    printf("Starting input processing.\n");
    thread_queue_run(tqueue);

    if (summary_stats_file)
        print_pair_stats(summary_stats_file);

    if (dist_fh) fclose(dist_fh);
    if (comp_fh) fclose(comp_fh);
    if (indel_fh) fclose(indel_fh);

    locus_diff_free(tqueue);

    if (opts.bf_par.readgroup_include_hash)
        free_readgroup_hash(opts.bf_par.readgroup_include_hash);
    
    printf("Finished.\n");

    return 0;
}
