#include <vector>
#include <utility>
#include <numeric>

#include "henry/tools.h"
#include "henry/error_estimate.h"
#include "henry/integrands.h"
#include "henry/sampling.h"
#include "henry/pileup_tools.h"
#include "henry/stats_tools.h"
#include "henry/stats_tools.h"

#include "alglib/lbfgs.h"

int main(int argc, char ** argv)
{

    char * base_qual_prior_file = argv[1];
    int min_quality_score = atoi(argv[2]);
    char * pileup_input_file = argv[3];
    char * posterior_output_file = argv[4];


    printf ("base_qual_prior_file: %s\n"
            "min_quality_score: %i\n"
            "pileup_input_file: %s\n"
            "posterior_output_file: %s\n",

            base_qual_prior_file,
            min_quality_score,
            pileup_input_file,
            posterior_output_file
            );
    

    FILE * pileup_input_fh = fopen(pileup_input_file, "r");
    if (pileup_input_fh == NULL)
    {
        fprintf(stderr, "Couldn't open pileup input file %s\n",
                pileup_input_file);
        exit(1);
    }

    FILE * posterior_output_fh = fopen(posterior_output_file, "w");
    if (posterior_output_fh == NULL)
    {
        fprintf(stderr, "Couldn't open posterior_output_file %s\n",
                posterior_output_file);
    }


    //we are integrating the actual posterior
    PileupSummary::initialize();

    char line_label[1000];
    size_t full_ndim = 4;


    PileupSummary summary(0);

    double mode_tolerance = 1e-30;
    size_t max_modefinding_iterations = 30;

    bool may_underflow = true;

    
    while (! feof(pileup_input_fh))
    {

        int next_char = fgetc(pileup_input_fh);
        if (next_char == EOF)
        {
            break;
        }
        ungetc(next_char, pileup_input_fh);

        bool succeeded = summary.load_line(pileup_input_fh);

        if (! succeeded)
        {
            fprintf(stderr, "Couldn't parse pileup line\n");
            exit(1);
        }
        
        PileupSummary const& p = summary;

        std::vector<BaseComp> base_observations =
            ExpandTo4State(p._bases_upper, p._quality_codes, p._read_depth, min_quality_score);
        
        BASE_COMP_COUNTS base_comp_counts = TallyBaseCompCounts(base_observations);

        size_t effective_depth = base_observations.size();

        ErrorEstimate error_estimate(base_comp_counts);
        error_estimate.Initialize(base_qual_prior_file);


        //the posterior may be called with 3 or 4 dimensions but it
        //ignores the 4th, assuming normalized coordinates
        Posterior posterior(&error_estimate, may_underflow, full_ndim);

        posterior.initialize(mode_tolerance, max_modefinding_iterations, 0);

        sprintf(line_label, "%s\t%i\t%c\t%i\t%Zu", summary._reference, 
                summary._position, summary._reference_base, summary._read_depth,
                effective_depth);
        
        for (size_t d = 0; d != 4; ++d)
        {
            fprintf(posterior_output_fh, "%s\t%c\t%10.8g\n", line_label,
                    ErrorEstimate::nucleotides[d], posterior.mode_point[d]);
        }
        fflush(posterior_output_fh);
    }
    
    fclose(pileup_input_fh);
    fclose(posterior_output_fh);

    return 0;

}
