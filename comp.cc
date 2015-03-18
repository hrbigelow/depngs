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
            "\nUsage: dep comp [options] input.jpd input.pileup contig_order.rdb output.comp [output.points]\n" 
            "Options:\n\n"
            "-l STRING   %s\n"
            "-r STRING   range file (lines 'contig<tab>start<tab>end') to process [blank] (blank = process whole file)\n"
            "-T INT      number of sample points used for tuning proposal distribution [1000]\n"
            "-t INT      number of threads to use [1]\n"
            "-m INT      number bytes of memory to use [%Zu]\n"
            "-f INT      number of sample points used for final quantiles estimation [1000]\n"
            "-y FLOAT    min mutation distance (0-1 scale) from ref to call as changed locus [0.2]\n"
            "-X FLOAT    confidence for -y [0.99]\n"
            "-Z FLOAT    confidence for binomial estimation [same as -X]\n"
            // "-a INT      target autocorrelation offset.  once reached, proposal tuning halts [6]\n"
            // "-I INT      maximum number of proposal tuning iterations [10]\n"
            // "-M FLOAT    autocorrelation maximum value [0.2]\n"
            "-C STRING   quantiles file [\"0.005 0.05 0.5 0.95 0.995\\n\"]\n"
            "-p FLOAT    alpha value for dirichlet prior [0.1]\n"
            "-q INT      %s\n"
            "-F STRING   Fastq offset type if known (one of Sanger, Solexa, Illumina13, Illumina15) [None]\n"
            "-v <empty>  if present, be verbose [absent]\n"
            "\n"
            "Output fields are:\n"
            "label_string algorithm reference position reference_base "
            "read_depth effective_depth full_anomaly_score pos_strand_anomaly_score "
            "neg_strand_anomaly_score inferred_base rank_order mean mode [quantiles]\n",

            Usage::label_string,
            1024 * 1024 * 1024 * 4l, // 4 GB memory
            Usage::quality_string
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
