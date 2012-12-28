/*
  compile with:
  make -C ~/arachne -B -k DEBUGGING=yes NO_DEPEND=yes PileupToPoints SYS_OPT=-I${HOME}/usr/include LINK_LIBS="-L${HOME}/usr/lib -lgsl -lgslcblas"
*/


#include <cstdlib>

#include "henry/tools.h"
#include "henry/pileup_tools.h"
#include "henry/error_estimate.h"




int main(int argc, char ** argv)
{

    char * pileup_input_file = argv[1];
    char * points_output_file = argv[2];

    //lookup table to convert a pileup code to reduced code for
    //initial processing
    
    FILE * points_fh = fopen(points_output_file, "w");
    if (points_fh == NULL)
    {
        fprintf(stderr, "Couldn't open points curve output file %s\n",
                points_output_file);
        exit(1);
    }

    FILE * pileup_input_fh = fopen(pileup_input_file, "r");
    if (pileup_input_fh == NULL)
    {
        fprintf(stderr, "Couldn't open pileup input file %s\n",
                pileup_input_file);
        exit(1);
    }

    PileupSummary::initialize();

    fprintf(points_fh, "locus\tA\tC\tG\tT\t\n");

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

        for (size_t b = 0; b != base_observations.size(); ++b)
        {
            BaseComp & bc = base_observations[b];
            fprintf(points_fh, "%i\t%f\t%f\t%f\t%f\n", p._position, bc.A, bc.C, bc.G, bc.T);
        }
    }

    fclose(pileup_input_fh);
    fclose(points_fh);

    return 0;
}
