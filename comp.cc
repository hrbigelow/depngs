#include <fstream>
#include <string.h>
#include <unistd.h>
#include <gsl/gsl_randist.h>

#include "run_comp.h"
#include "comp_worker.h"

#include "usage_strings.h"

int comp_usage()
{
    fprintf(stderr, 
            "\nUsage: dep comp [options] samples.rdb contig_order.rdb\n" 
            "Options:\n\n"
            "FILES\n"
            "-d STRING   (optional) name of output distance file.  If absent, do not perform distance calculation.\n"
            "-c STRING   (optional) name of output composition file.  If absent, do not output composition marginal estimates.\n"
            "-i STRING   (optional) name of output indel composition file.  If absent, do not perform indel distance calculation.\n"
            "-r STRING   input range file (lines 'contig<tab>start<tab>end') to process [blank] (blank = process whole file)\n"
            "-x STRING   (optional) name of output summary statistics file of differing loci.\n"
            "\n"
            "STATISTICAL PARAMETERS\n"
            "-y FLOAT    MIN_DIST, min mutation distance (0-1 scale) to call as changed [0.2]\n"
            "-X FLOAT    POST_CONF, confidence for -y [0.99]\n"
            "-Z FLOAT    BETA_CONF, confidence for binomial estimation [0.9999]\n"
            "-f INT      MAX_POINTS, max # of sample points for binomial test [10000]\n"
            "-p FLOAT    PRIOR_ALPHA, alpha value for dirichlet prior [0.1]\n"
            "\n"
            "OTHER\n"
            "-q INT      minimum quality score to include bases as evidence [5]\n"
            "-F STRING   Fastq offset type if known (one of Sanger, Solexa, Illumina13, Illumina15) [None]\n"
            "-g <empty>  if present, print extra pileup fields in the output distance file [absent]\n"
            "-t INT      number of threads to use [1]\n"
            "-R INT      number of readers to use [2] (See NOTE)\n"
            "-m INT      number bytes of memory to use [50e9]\n"
            "\n"
            "These options affect how sampling is done.\n"
            "\n"
            "-Q STRING   quantiles string [\"0.005,0.05,0.5,0.95,0.995\"]\n"
            "\n"
            "samples.rdb has lines of <sample_id><tab></path/to/sample.jpd><tab></path/to/sample.pileup>\n"
            "sample_pairings.rdb has lines of <sample_id><tab><sample_id>\n"
            "defining which pairs of samples are to be compared.  <sample_id> correspond\n"
            "with those in samples.rdb\n"
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



int main_comp(int argc, char ** argv)
{

    char label_string[100];
    strcpy(label_string, "comp");

    size_t num_threads = 1;
    size_t max_mem = 1024l * 1024l * 1024l * 4l;
    float prior_alpha = 0.1;

    struct posterior_settings pset = {
        { prior_alpha, prior_alpha, prior_alpha, prior_alpha },
        1000,    /* max_sample_points */
        0.2,     /* min_dist */
        0.99,    /* post_confidence */
        0.99,    /* beta_confidence */
        5       // min_quality_score
    };


    /* initialize logu with suitably distributed uniform values.  this
       is cut-and-pasted from dist.cc.  Likely this logu scheme will
       be replaced by a non-MH sampling scheme.
     */
    // pset.logu = (double *)malloc(sizeof(double) * pset.final_n_points);
    // double n_inv = 1.0 / (double)(pset.final_n_points);

    // unsigned p;
    // for (p = 0; p != pset.final_n_points; ++p)
    //     pset.logu[p] = log(p * n_inv);

    // gsl_rng *randgen = gsl_rng_alloc(gsl_rng_taus);
    // gsl_ran_shuffle(randgen, pset.logu, pset.final_n_points, sizeof(pset.logu[0]));
    // gsl_rng_free(randgen);


    // double test_quantile = 0.01;
    // double min_test_quantile_value = 0;

    char quantiles_file[100];
    strcpy(quantiles_file, "/dev/null");

    const char *fastq_type = NULL;

    const char 
        *jpd_data_params_file,
        *pileup_input_file,
        *contig_order_file,
        *query_range_file = NULL,
        *posterior_output_file,
        *cdfs_output_file;

    bool verbose = false;

    char c;
    while ((c = getopt(argc, argv, "l:r:t:m:f:y:X:Z:C:p:q:F:v")) >= 0)
    {
        switch(c)
        {
        case 'l': strcpy(label_string, optarg); break;
        case 'r': query_range_file = optarg; break;
        case 't': num_threads = (size_t)atof(optarg); break;
        case 'm': max_mem = (size_t)atof(optarg); break;
        // case 'T': pset.tuning_n_points = (size_t)atof(optarg); break;
        case 'f': pset.max_sample_points = (size_t)atof(optarg); break;
        case 'y': pset.min_dist = atof(optarg); break;
        case 'X': pset.post_confidence = atof(optarg); break;
        case 'Z': pset.beta_confidence = atof(optarg); break;
        // case 'a': pset.target_autocor_offset = (size_t)atof(optarg); break;
        // case 'I': pset.max_tuning_iterations = (size_t)atof(optarg); break;
        // case 'M': pset.autocor_max_value = atof(optarg); break;
        case 'C': strcpy(quantiles_file, optarg); break;
        case 'p': prior_alpha = atof(optarg); break;
        case 'q': pset.min_quality_score = static_cast<size_t>(atoi(optarg)); break;
        case 'F': fastq_type = optarg; break;
        case 'v': verbose = true; break;
        default: return comp_usage(); break;
        }
    }
    if (argc - optind < 4 || argc - optind > 5)
        return comp_usage();

    jpd_data_params_file = argv[optind];
    pileup_input_file = argv[optind + 1];
    contig_order_file = argv[optind + 2];
    posterior_output_file = argv[optind + 3];

    cdfs_output_file = (optind + 4 < argc) ? argv[optind + 4] : "/dev/null";

    int rval = run_comp(max_mem,
                        num_threads,
                        fastq_type,
                        label_string,
                        quantiles_file,
                        pileup_input_file,
                        contig_order_file,
                        query_range_file,
                        jpd_data_params_file,
                        posterior_output_file,
                        cdfs_output_file,
                        &pset,
                        verbose);

    return rval;
}
