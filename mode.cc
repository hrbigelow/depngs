#include <vector>
#include <utility>
#include <numeric>
#include <unistd.h>

#include <gsl/gsl_math.h>
#include <gsl/gsl_statistics_double.h>

#include <sys/timeb.h>

#include "tools.h"
#include "error_estimate.h"
#include "integrands.h"
#include "metropolis.h"
#include "sampling.h"
#include "pileup_tools.h"
#include "stats_tools.h"
#include "dirichlet.h"
#include "slice_sampling.h"
#include "stats_tools.h"
#include "nucleotide_stats.h"
#include "simulation.h"
#include "anomaly_tools.h"
#include "transformation.h"
#include "samutil/file_utils.h"

#include "usage_strings.h"

int mode_usage()
{
    fprintf(stderr, 
            "\nUsage: dep mode [options] input.jpd input.pileup output.mode\n" 
            "Options:\n\n"
            "-l STRING   %s\n"
            "-p FILE     alpha values for dirichlet prior [\"0.1 0.1 0.1 0.1\\n\"]\n"
            "-q INT      %s\n"
            "-v <empty>  if present, be verbose [absent]\n"
            "-a <empty>  if present, compute anomaly scores [absent]\n"
            "\n"
            "Output fields are:\n"
            "label_string reference position reference_base "
            "read_depth effective_depth full_anomaly_score pos_strand_anomaly_score "
            "neg_strand_anomaly_score modeA modeC modeG modeT\n",
            Usage::mode_label_string, Usage::quality_string
            );
    return 1;
}

extern char *optarg;
extern int optind;



int main_mode(int argc, char ** argv)
{

    char label_string[100];
    strcpy(label_string, "mode");

    char prior_alphas_file[100];
    strcpy(prior_alphas_file, "/dev/null");

    size_t min_quality_score = 5;

    double default_prior_alpha = 0.1;

    char const* jpd_data_params_file;
    char const* pileup_input_file;
    char const* posterior_output_file;

    bool verbose = false;
    bool compute_anomaly = false;
    size_t max_mem = 1024 * 1024 * 1024;

    char c;
    while ((c = getopt(argc, argv, "d:l:p:q:vam:")) >= 0)
    {
        switch(c)
        {
        case 'l': strcpy(label_string, optarg); break;
        case 'p': strcpy(prior_alphas_file, optarg); break;
        case 'q': min_quality_score = static_cast<size_t>(atoi(optarg)); break;
        case 'v': verbose = true; break;
        case 'a': compute_anomaly = true; break;
        case 'm': max_mem = static_cast<size_t>(atoi(optarg)); break;
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

    double * prior_alphas;
    size_t num_prior_alphas;

    if (strcmp(prior_alphas_file, "/dev/null") == 0)
    {
        prior_alphas = new double[4];
        std::fill(prior_alphas, prior_alphas + 4, default_prior_alpha);
    }
    else
    {
        prior_alphas = ParseNumbersFile(prior_alphas_file, & num_prior_alphas);
    }

    gsl_rng * rand_gen = gsl_rng_alloc(gsl_rng_taus);
    timeb millitime;
    ftime(& millitime);
    gsl_rng_set(rand_gen, millitime.millitm);
    
    size_t full_ndim = 4;

    FILE * posterior_output_fh = fopen(posterior_output_file, "w");
    if (posterior_output_fh == NULL)
    {
        fprintf(stderr, "Couldn't open posterior_output_file %s\n",
                posterior_output_file);
    }

    FILE * pileup_input_fh = open_if_present(pileup_input_file, "r");

    PileupSummary pileup(0);
    size_t chunk_size = max_mem;
    char * chunk_buffer_in = new char[chunk_size + 1];

    FastqType ftype = pileup.FastqFileType(pileup_input_file, chunk_buffer_in, chunk_size);

    if (ftype == None)
    {
        fprintf(stderr, "Error: Couldn't determine quality scale for pileup input file %s\n",
                pileup_input_file);
        exit(1);
    }

    PileupSummary::SetFtype(ftype);
   

    //we are integrating the actual posterior
    double mode_tolerance = 1e-60;
    size_t max_modefinding_iterations = 3000;
    double initial_point[] = { 0.25, 0.25, 0.25, 0.25 };


    bool may_underflow = true;

    size_t effective_depth;
    
    ErrorEstimate model;
    model.set_composition_prior_alphas(prior_alphas);
    
    //1. initialize NucleotideStats from JPD_DATA
    NucleotideStats params;
    params.initialize(jpd_data_params_file);

    //3. Set model parameters
    model.model_params = & params;

    //4. Construct posterior
    Posterior posterior(&model, may_underflow, full_ndim);

    size_t nbytes_read, nbytes_unused = 0;
    char * last_fragment;
    char * read_pointer = chunk_buffer_in;
    
    while (! feof(pileup_input_fh))
    {
        nbytes_read = fread(read_pointer, 1, chunk_size - nbytes_unused, pileup_input_fh);

        std::vector<char *> pileup_lines =
            FileUtils::find_complete_lines_nullify(chunk_buffer_in, & last_fragment);

        std::vector<char *>::iterator pit;
        read_pointer[nbytes_read] = '\0';

        for (pit = pileup_lines.begin(); pit != pileup_lines.end(); ++pit)
        {

            PileupSummary locus(0);
            locus.load_line((*pit));
            locus.parse(min_quality_score);
            params.pack(& locus.counts);
            
            //divide locus data to plus and minus-strand data
            
            effective_depth = locus.read_depth;

            posterior.model()->locus_data = & locus.counts;

            posterior.initialize(mode_tolerance, max_modefinding_iterations, 
                                 initial_point, verbose);

        
            fprintf(posterior_output_fh, 
                    "%s\t%s\t%i\t%c\t%Zu\t%Zu\t%5.5f\t%5.5f\t%5.5f\t%5.5f",
                    label_string, 
                    locus.reference, 
                    locus.position, 
                    locus.reference_base, 
                    locus.read_depth,
                    effective_depth,
                    posterior.mode_point[0],
                    posterior.mode_point[1],
                    posterior.mode_point[2],
                    posterior.mode_point[3]
                    );

            if (compute_anomaly)
            {
                // double pos_anomaly_score =
                //     strand_locus_anomaly_score(posterior, global_counts,
                //                                locus, data_reader, '+', verbose);
                
                // double neg_anomaly_score =
                //     strand_locus_anomaly_score(posterior, global_counts,
                //                                locus, data_reader, '-', verbose);
                double full_anomaly_score =
                    relative_entropy_anomaly(params.cpd_buffer, & locus.counts, posterior.mode_point);

                // fprintf(posterior_output_fh,
                //         "\t%5.5lf\t%5.5lf\t%5.5lf\n",
                //         full_anomaly_score,
                //         pos_anomaly_score,
                //         neg_anomaly_score);
            }
            else
            {
                fprintf(posterior_output_fh, "\n");
            }
                        
            //fflush(posterior_output_fh);

        }
        nbytes_unused = strlen(last_fragment);
        memmove(chunk_buffer_in, last_fragment, nbytes_unused);
        read_pointer = chunk_buffer_in + nbytes_unused;
    }

    fclose(pileup_input_fh);
    delete chunk_buffer_in;
    
    fclose(posterior_output_fh);
    
    gsl_rng_free(rand_gen);
    
    delete prior_alphas;
    return 0;
}
