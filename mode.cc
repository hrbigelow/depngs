#include <string.h>

#include "run_comp_or_mode.h"
#include "comp_functor.h"

#include "usage_strings.h"

int mode_usage()
{
    fprintf(stderr, 
            "\nUsage: dep mode [options] input.jpd input.pileup output.mode\n" 
            "Options:\n\n"
            "-l STRING   %s\n"
            "-t INT      number of threads to use [1]\n"
            "-z REAL     gradient tolerance used in finding the posterior mode (gsl_multimin) [1e-5]\n"
            "-m INT      number bytes of memory to use [%Zu]\n"
            "-p FILE     alpha values for dirichlet prior [\"0.1 0.1 0.1 0.1\\n\"]\n"
            "-q INT      %s\n"
            "-T INT      number of sample points used for tuning proposal distribution [1000]\n"
            "-f INT      number of sample points used for final quantiles estimation [10000]\n"
            "-v <empty>  if present, be verbose [absent]\n"
            "\n"
            "Output fields are:\n"
            "label_string reference position reference_base "
            "read_depth effective_depth full_anomaly_score pos_strand_anomaly_score "
            "neg_strand_anomaly_score modeA modeC modeG modeT\n",
            Usage::mode_label_string, 
            1024 * 1024 * 1024 * 4l, // 4 GB memory
            Usage::quality_string
            );
    return 1;
}

extern char *optarg;
extern int optind;



int main_mode(int argc, char ** argv)
{

    char label_string[100];
    strcpy(label_string, "mode");

    size_t num_threads = 1;
    size_t max_mem = 1024l * 1024l * 1024l * 4l;
    double gradient_tolerance = 1e-5;

    char prior_alphas_file[100];
    strcpy(prior_alphas_file, "/dev/null");

    size_t min_quality_score = 5;

    char const* jpd_data_params_file;
    char const* pileup_input_file;
    char const* posterior_output_file;

    size_t tuning_num_points = 1e3;
    size_t final_num_points = 1e4;

    bool verbose = false;

    char c;
    while ((c = getopt(argc, argv, "l:t:m:z:p:qT:f::v")) >= 0)
    {
        switch(c)
        {
        case 'l': strcpy(label_string, optarg); break;
        case 't': num_threads = static_cast<size_t>(atof(optarg)); break;
        case 'm': max_mem = static_cast<size_t>(atof(optarg)); break;
        case 'z': gradient_tolerance = atof(optarg); break;
        case 'p': strcpy(prior_alphas_file, optarg); break;
        case 'q': min_quality_score = static_cast<size_t>(atoi(optarg)); break;
        case 'T': tuning_num_points = static_cast<size_t>(atof(optarg)); break;
        case 'f': final_num_points = static_cast<size_t>(atof(optarg)); break;
        case 'v': verbose = true; break;
        default: return mode_usage(); break;
        }
    }
    if (argc - optind != 3)
    {
        return mode_usage();
    }

    jpd_data_params_file = argv[optind];
    pileup_input_file = argv[optind + 1];
    posterior_output_file = argv[optind + 2];

    return run_comp_or_mode(max_mem,
                            num_threads,
                            min_quality_score,
                            label_string,
                            "/dev/null",
                            prior_alphas_file,
                            pileup_input_file,
                            jpd_data_params_file,
                            posterior_output_file,
                            "/dev/null",
                            gradient_tolerance,
                            tuning_num_points,
                            final_num_points,
                            verbose,
                            &mode_worker);

}
