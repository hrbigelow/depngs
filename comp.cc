#include <fstream>
#include <string.h>

#include "run_comp_or_mode.h"
#include "comp_functor.h"

#include "usage_strings.h"

int comp_usage()
{
    fprintf(stderr, 
            "\nUsage: dep comp [options] input.jpd input.pileup output.comp [output.points]\n" 
            "Options:\n\n"
            "-l STRING   %s\n"
            "-T INT      number of sample points used for tuning proposal distribution [1000]\n"
            "-t INT      number of threads to use [1]\n"
            "-m INT      number bytes of memory to use [%Zu]\n"
            "-z REAL     gradient tolerance used in finding the posterior mode (gsl_multimin) [1e-5]\n"
            "-f INT      number of sample points used for final quantiles estimation [1000]\n"
            "-X FLOAT    test quantile used for -y [0.01]\n"
            "-y FLOAT    minimum test quantile value needed to output locus composition [0]\n"
            "-a INT      target autocorrelation offset.  once reached, proposal tuning halts [6]\n"
            "-i INT      maximum number of proposal tuning iterations [10]\n"
            "-s INT      number of loci simulations for estimating anomaly score (if zero, no estimate provided) [0]\n"
            "-M FLOAT    autocorrelation maximum value [6]\n"
            "-C STRING   quantiles file [\"0.005 0.05 0.5 0.95 0.995\\n\"]\n"
            "-p STRING   alpha values for dirichlet prior [\"0.1 0.1 0.1 0.1\\n\"]\n"
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

    size_t tuning_num_points = 1e3;
    size_t final_num_points = 1e3;

    double test_quantile = 0.01;
    double min_test_quantile_value = 0;

    size_t target_autocor_offset = 6;
    size_t max_tuning_iterations = 10;
    double gradient_tolerance = 1e-5;

    double autocor_max_value = 6;



    size_t nsim_loci = 0;

    char quantiles_file[100];
    strcpy(quantiles_file, "/dev/null");

    char prior_alphas_file[100];
    strcpy(prior_alphas_file, "/dev/null");

    size_t min_quality_score = 5;
    const char *fastq_type = "None";

    char const* jpd_data_params_file;
    char const* pileup_input_file;
    char const* posterior_output_file;
    char const* cdfs_output_file;

    bool verbose = false;

    char c;
    while ((c = getopt(argc, argv, "l:t:m:z:T:f:X:y:a:s:i:M:C:p:q:F:v")) >= 0)
    {
        switch(c)
        {
        case 'l': strcpy(label_string, optarg); break;
        case 't': num_threads = static_cast<size_t>(atof(optarg)); break;
        case 'm': max_mem = static_cast<size_t>(atof(optarg)); break;
        case 'z': gradient_tolerance = atof(optarg); break;
        case 'T': tuning_num_points = static_cast<size_t>(atof(optarg)); break;
        case 'f': final_num_points = static_cast<size_t>(atof(optarg)); break;
        case 'X': test_quantile = atof(optarg); break;
        case 'y': min_test_quantile_value = atof(optarg); break;
        case 'a': target_autocor_offset = static_cast<size_t>(atof(optarg)); break;
        case 'i': max_tuning_iterations = static_cast<size_t>(atof(optarg)); break;
        case 's': nsim_loci = static_cast<size_t>(atof(optarg)); break;
        case 'M': autocor_max_value = atof(optarg); break;
        case 'C': strcpy(quantiles_file, optarg); break;
        case 'p': strcpy(prior_alphas_file, optarg); break;
        case 'q': min_quality_score = static_cast<size_t>(atoi(optarg)); break;
        case 'F': fastq_type = optarg; break;
        case 'v': verbose = true; break;
        default: return comp_usage(); break;
        }
    }
    if (argc - optind != 3)
    {
        return comp_usage();
    }

    jpd_data_params_file = argv[optind];
    pileup_input_file = argv[optind + 1];
    posterior_output_file = argv[optind + 2];

    cdfs_output_file = (optind + 3 < argc) ? argv[optind + 3] : "/dev/null";

    return run_comp_or_mode(max_mem,
                            num_threads,
                            min_quality_score,
                            fastq_type,
                            label_string,
                            quantiles_file,
                            prior_alphas_file,
                            pileup_input_file,
                            jpd_data_params_file,
                            posterior_output_file,
                            cdfs_output_file,
                            gradient_tolerance,
                            tuning_num_points,
                            final_num_points,
                            test_quantile,
                            min_test_quantile_value,
                            verbose,
                            &comp_worker);

}
