/*
  compile with:

  make -C ~/broadcrd -B -k DEBUGGING=yes NO_DEPEND=yes ModeSampleComposition SYS_OPT=-I${HOME}/usr/include LINK_LIBS="-L${HOME}/usr/lib -lgsl -lgslcblas -lm -lgmpxx -lgmp"

  make -C ~/broadcrd -B -k NO_DEPEND=yes ModeSampleComposition SYS_OPT="-O3 -I${HOME}/usr/include" LINK_LIBS="-L${HOME}/usr/lib -lgsl -lgslcblas -lm -lgmpxx -lgmp"
*/

#include <cstdlib>
#include <cstdio>
#include <set>
#include <string>
#include <cstring>
#include <cmath>
#include <cassert>
#include <algorithm>
#include <inttypes.h>

#include <gsl/gsl_vector.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_multimin.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_monte.h>
#include <gsl/gsl_monte_plain.h>
#include <gsl/gsl_monte_miser.h>
#include <gsl/gsl_monte_vegas.h>

#include "henry/error_estimate.h"
#include "henry/tools.h"
#include "henry/pileup_tools.h"
#include "henry/integrands.h"




int main(int argc, char ** argv)
{

    double mode_tolerance = atof(argv[1]);
    double scale_increment = atof(argv[2]);
    int min_quality_score = atoi(argv[3]);

    fprintf(stderr, 
            "mode_tolerance: %g\n"
            "scale_increment: %g\n"
            "min_quality_score: %i\n",
            mode_tolerance,
            scale_increment,
            min_quality_score);


   
    PileupSummary::initialize();

    while (! feof(stdin))
    {

        int next_char = fgetc(stdin);
        if (next_char == EOF)
        {
            break;
        }
        ungetc(next_char, stdin);

        PileupSummary summary(0);

        summary.load_line(stdin);
        
        PileupSummary const& p = summary;

        std::vector<BaseComp> base_observations =
            ExpandTo4State(p._bases_upper, p._quality_codes, p._read_depth, min_quality_score);

        BASE_COMP_COUNTS base_comp_counts = TallyBaseCompCounts(base_observations);

        ErrorEstimate error_estimate(base_comp_counts);

        //Find the mode point
        bool is_log2_form = false;
        size_t num_dimensions = 3;

        Posterior posterior(&error_estimate, is_log2_form, num_dimensions);

        size_t max_function_evals = 100000;

        posterior.initialize(mode_tolerance, max_function_evals, 0);

        int const* bi = PileupSummary::base_to_index;
        int const* bc = summary.base_counts;

        for (size_t base = 0; base != 4; ++base)
        {

            printf("%s\t%i\t%i\t%i\t%i\t%i\t%i\t",
                   summary._reference,
                   summary._position, summary._read_depth,
                   bc[bi['A']] + bc[bi['a']],
                   bc[bi['C']] + bc[bi['c']],
                   bc[bi['G']] + bc[bi['g']],
                   bc[bi['T']] + bc[bi['t']]);
            
            printf("%c\t%g\n", 
                    ErrorEstimate::nucleotides[base],
                    posterior.mode_point[base]);
        }

        fflush(stdout);
    }
    
    //fclose(stdin);
    
    return 0;
}
