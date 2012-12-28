/*
  compile with:

  make -C ~/arachne -B -k DEBUGGING=yes NO_DEPEND=yes GridSampleComposition SYS_OPT=-I${HOME}/usr/include LINK_LIBS="-L${HOME}/usr/lib -lgsl -lgslcblas -lm"

  make -C ~/arachne -B -k NO_DEPEND=yes GridSampleComposition SYS_OPT="-O3 -I${HOME}/usr/include" LINK_LIBS="-L${HOME}/usr/lib -lgsl -lgslcblas -lcuba -lm"
*/

#include <cstdlib>
#include <cstdio>
#include <set>
#include <string>
#include <cmath>
#include <cassert>

#include <gsl/gsl_vector.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_multimin.h>
#include <gsl/gsl_math.h>

#include "henry/error_estimate.h"
#include "henry/tools.h"
#include "henry/pileup_tools.h"





int main(int argc, char ** argv)
{

    double mode_tolerance = atof(argv[1]);
    int npoints = atoi(argv[2]);
    double grid_spacing = atof(argv[3]);
    char * points_output_file = argv[4];
    char * grid_output_file = argv[5];
    char * pileup_input_file = argv[6];

    printf ("mode_tolerance: %g\n"
            "npoints: %i\n"
            "grid_spacing: %g\n"
            "points_output_file: %s\n"
            "grid_output_file: %s\n"
            "pileup_input_file: %s\n",
            mode_tolerance,
            npoints,
            grid_spacing,
            points_output_file,
            grid_output_file,
            pileup_input_file);


    //lookup table to convert a pileup code to reduced code for
    //initial processing
    
    FILE * points_fh = fopen(points_output_file, "w");
    if (points_fh == NULL)
    {
        fprintf(stderr, "Couldn't open points curve output file %s\n",
                points_output_file);
        exit(1);
    }

    FILE * grid_fh = fopen(grid_output_file, "w");
    if (grid_fh == NULL)
    {
        fprintf(stderr, "Couldn't open grid curve output file %s\n",
                grid_output_file);
        exit(1);
    }

    FILE * pileup_input_fh = fopen(pileup_input_file, "r");
    if (pileup_input_fh == NULL)
    {
        fprintf(stderr, "Couldn't open pileup input file %s\n",
                pileup_input_file);
        exit(1);
    }

    int last_index;
    int scanned;

    ErrorEstimate error_estimate;
   
    int scanned_fields;

    PileupSummary::initialize();

    while (! feof(pileup_input_fh))
    {

        int next_char = fgetc(pileup_input_fh);
        if (next_char == EOF)
        {
            break;
        }
        ungetc(next_char, pileup_input_fh);

        PileupSummary summary(0);

        summary.load_line(pileup_input_fh);
        
        PileupSummary const& p = summary;

        std::vector<BaseComp> base_observations =
            ExpandTo4State(p._bases_upper, p._quality_codes, p._read_depth);

        BASE_COMP_COUNTS base_comp_counts = TallyBaseCompCounts(base_observations);

        double mode[3];

        size_t num_required_iterations =
            error_estimate.find_mode_point(mode_tolerance, 
                                           max_iterations,
                                           mode);
        
        BaseComp mode_comp(mode[0], mode[1], mode[2], 1.0 - mode[0] - mode[1] - mode[2]);
        
        //deallocate minimizer and x vector
        gsl_multimin_fminimizer_free(minimizer);
        gsl_vector_free(x);
        gsl_vector_free(ss);

        REAL log2_mode = error_estimate.Log2Posterior(base_comp_counts, mode_comp);

        //print posterior values for a cartesian grid about the mode, using
        //nx, ny, nz points and a given grid spacing

        REAL xlow[3];
        for (int dim = 0; dim != 3; ++dim)
        {
            xlow[dim] = mode[dim] - (npoints / 2) * grid_spacing;
        }

        for (int x = 0; x != npoints; ++x)
        {
            REAL xc = xlow[0] + x * grid_spacing;
            for (int y = 0; y != npoints; ++y)
            {
                REAL yc = xlow[1] + y * grid_spacing;
                for (int z = 0; z != npoints; ++z)
                {
                    REAL zc = xlow[2] + z * grid_spacing;
                    BaseComp comp(xc, yc, zc, 1.0 - xc - yc - zc);
                    REAL fval = 
                        error_estimate.ScaledPosterior(base_comp_counts, comp, log2_mode);

                    REAL log_fval =
                        error_estimate.Log2Posterior(base_comp_counts, comp);

                    fprintf(points_fh, "%20.20f\t%20.20f\t%20.20f\t%20.20f\t%20.20g\n",
                            xc, yc, zc, fval, log_fval - log2_mode);
                }
            }
        }


        for (int x = 0; x != npoints; ++x)
        {
            REAL xc = xlow[0] + x * grid_spacing;
            REAL yc = xlow[1] + x * grid_spacing;
            REAL zc = xlow[2] + x * grid_spacing;
            fprintf(grid_fh, "%20.20f\t%20.20f\t%20.20f\n", xc, yc, zc);
        }
    }

    fclose(pileup_input_fh);
    fclose(points_fh);
    fclose(grid_fh);

    return 0;
}
