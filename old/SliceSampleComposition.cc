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
#include <gsl/gsl_math.h>


#include "henry/error_estimate.h"
#include "henry/tools.h"
#include "henry/pileup_tools.h"
#include "henry/slice_sampling.h"
#include "henry/integrands.h"




int main(int argc, char ** argv)
{

    char * mass_fractions_file = argv[1];
    char * posterior_output_file = argv[2];
    char * points_output_file = argv[3];
    char * pileup_input_file = argv[4];
    size_t number_sample_points = static_cast<size_t>(atof(argv[5]));
    size_t initial_sampling_range = static_cast<size_t>(atof(argv[6]));
    int print_points = atoi(argv[7]);
    int min_quality_score = atoi(argv[8]);

    printf ("mass_fractions_file: %s\n"
            "posterior_output_file: %s\n"
            "points_output_file: %s\n"
            "pileup_input_file: %s\n"
            "number_sample_points: %i\n"
            "initial_sampling_range: %i\n"
            "print_points: %i\n"
            "min_quality_score: %i\n",

            mass_fractions_file,
            posterior_output_file,
            points_output_file,
            pileup_input_file,
            static_cast<int>(number_sample_points),
            static_cast<int>(initial_sampling_range),
            print_points,
            min_quality_score);

    int num_dimensions = 3; //nucleotide space
    int num_bits_per_dim = 62;
    bool is_log2_integrand = true;

    float mass_fractions[100];
    FILE * mass_fractions_fh = fopen(mass_fractions_file, "r");
    if (mass_fractions_fh == NULL)
    {
        fprintf(stderr, "Couldn't open integration points file %s\n", 
                mass_fractions_file);
        exit(1);
    }

    int b = 0;
    while (! feof(mass_fractions_fh))
    {
        fscanf(mass_fractions_fh, "%f ", &mass_fractions[b]);
        ++b;
    }
    size_t num_mass_fractions = b;
    fclose(mass_fractions_fh);


    FILE * posterior_fh = fopen(posterior_output_file, "w");
    if (posterior_fh == NULL)
    {
        fprintf(stderr, "Couldn't open posterior curve output file %s\n"
                "Will not print out posterior curve sample points.\n",
                posterior_output_file);
    }

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


    //PileupSummary::initialize();
    SliceSampling slice_sampling(num_dimensions, num_bits_per_dim,
                                 is_log2_integrand);
    
    slice_sampling.Initialize();

    bool is_log2_form = false;
    double mode_tolerance = 1e-10;

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
        //int sample_alphabet_size = strlen(sample_base_alphabet);

        std::vector<BaseComp> base_observations =
            ExpandTo4State(p._bases, p._quality_codes, p._read_depth, min_quality_score);

        BASE_COMP_COUNTS base_comp_counts = TallyBaseCompCounts(base_observations);

        ErrorEstimate error_estimate(base_comp_counts);

        Posterior integrand(error_estimate, is_log2_form, num_dimensions);

        integrand.initialize(mode_tolerance, number_sample_points);

        //use slice sampling to get a set of sample points

        WEIGHTED_SAMPLE_MAP sample_points =
            slice_sampling.sample(&integrand, 
                                  integrand.mode_point, 
                                  initial_sampling_range, 
                                  number_sample_points);

        WEIGHTED_SAMPLE_MAP::iterator sample_iter;

        std::vector<WeightedSample *> sample_points_vec;

        for (sample_iter = sample_points.begin();
             sample_iter != sample_points.end();
             ++sample_iter)
        {
            sample_points_vec.push_back(&(*sample_iter).second);
        }
        
        size_t print_limit = 10000;

        if (print_limit > 0)
        if(0)
        {
            printf("\n********************************************************\n"
                   "%i samples\n"
                   "********************************************************\n",
                   static_cast<int>(sample_points.size()));
        }
        
        size_t sample_num = 0;
        for (size_t s = 0; s != sample_points_vec.size(); ++s)
        {

            if (sample_num >= print_limit)
            {
                break;
            }
            WeightedSample & w = *sample_points_vec[s];
            printf("%7i %14.14f %14.14f %14.14f %14.14f: val: %14.14f wt: %14.14f\n",
                   static_cast<int>(sample_num), w.x[0], w.x[1], w.x[2], w.x[3], w.val, w.weight);
            ++sample_num;
        }


        
        //find desired marginal integral boundaries
        std::vector<REAL> bound_sums(num_mass_fractions);
        std::fill(bound_sums.begin(), bound_sums.end(), 0.0);

        int const* bi = PileupSummary::base_to_index;
        int const* bc = summary.base_counts;

        for (size_t base = 0; base != 4; ++base)
        {

            fprintf(posterior_fh, "%i\t%i\t%i\t%i\t%i\t%i\t",
                    summary._position, summary._read_depth,
                    bc[bi['A']] + bc[bi['a']],
                    bc[bi['C']] + bc[bi['c']],
                    bc[bi['G']] + bc[bi['g']],
                    bc[bi['T']] + bc[bi['t']]);

            fprintf(posterior_fh, 
                    "%7i\t%c", static_cast<int>(sample_points.size()),
                    ErrorEstimate::nucleotides[base]);
                
            MarginalCumulativeDistribution(base, &sample_points_vec);

            for (size_t mi = 0; mi != num_mass_fractions; ++mi)
            {

                bool is_lower_bound = mass_fractions[mi] < 0.5;
                REAL bound = 
                    FindIntegralBound(&sample_points_vec, base,
                                      mass_fractions[mi],
                                      is_lower_bound);
                bound_sums[mi] += bound;

                fprintf(posterior_fh, "\t% 8.5f", bound);
            }
            fprintf(posterior_fh, "\n");
        }

        fprintf(posterior_fh, "%i\t%i\t%i\t%i\t%i\t%i\t",
                summary._position, summary._read_depth,
                bc[bi['A']] + bc[bi['a']],
                bc[bi['C']] + bc[bi['c']],
                bc[bi['G']] + bc[bi['g']],
                bc[bi['T']] + bc[bi['t']]);

        fprintf(posterior_fh, "%7i\t%c", 
                static_cast<int>(sample_points.size()),
                '+');

        for (size_t mi = 0; mi != num_mass_fractions; ++mi)
        {
            fprintf(posterior_fh, "\t% 8.5f", bound_sums[mi]);
        }
        fprintf(posterior_fh, "\n");


        fflush(posterior_fh);


        if (print_points)
        {
            fprintf(points_fh, 
                    "locus\tA\tC\tG\tT\tcA\tcC\tcG\tcT\t"
                    "integrand\tpoint_weight\n");
                    
            for (size_t s = 0; s != sample_points_vec.size(); ++s)
//             for (sample_iter = sample_points.begin();
//                  sample_iter != sample_points.end();
//                  ++sample_iter)
            {
                
                WeightedSample & w = *sample_points_vec[s];

                fprintf(points_fh, 
                        "%i\t%.10g\t%.10g\t%.10g\t%.10g"
                        "\t%.10g\t%.10g\t%.10g\t%.10g"
                        "\t%.10g\t%.10g\n",
                        summary._position, 
                        w.x[0], w.x[1], w.x[2], w.x[3], 
                        w.cdf[0], w.cdf[1], w.cdf[2], w.cdf[3],
                        w.val, w.weight);
            }
        }
    }
    
    fclose(pileup_input_fh);
    fclose(posterior_fh);
    fclose(points_fh);
    
    return 0;
}
