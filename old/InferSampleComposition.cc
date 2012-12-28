/*
  compile with:

  make -C ~/arachne -B -k DEBUGGING=yes NO_DEPEND=yes InferSampleComposition SYS_OPT=-I${HOME}/usr/include LINK_LIBS="-L${HOME}/usr/lib -lgsl -lgslcblas -lcuba_debug -lm -lgmpxx -lgmp"

  make -C ~/arachne -B -k NO_DEPEND=yes InferSampleComposition SYS_OPT="-O3 -I${HOME}/usr/include" LINK_LIBS="-L${HOME}/usr/lib -lgsl -lgslcblas -lcuba -lm -lgmpxx -lgmp"
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

#include <cuba.h>

#include <gsl/gsl_vector.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_multimin.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_monte.h>
#include <gsl/gsl_monte_plain.h>
#include <gsl/gsl_monte_miser.h>
#include <gsl/gsl_monte_vegas.h>

//#include "MainTools.h"

#include "henry/error_estimate.h"
#include "henry/tools.h"
#include "henry/pileup_tools.h"





int main(int argc, char ** argv)
{

    size_t max_function_evals = static_cast<size_t>(atof(argv[1]));
    size_t min_function_evals = static_cast<size_t>(atof(argv[2]));
    float epsilon_abs = atof(argv[3]);
    float epsilon_rel = atof(argv[4]);
    double mode_tolerance = atof(argv[5]);
    float fraction_used_points = atof(argv[6]);
    double scale_increment = atof(argv[7]);
    char * mass_fractions_file = argv[8];
    char * posterior_output_file = argv[9];
    char * points_output_file = argv[10];
    char * pileup_input_file = argv[11];
    int selected_alg = atoi(argv[12]);
    int collect_points = atoi(argv[13]);
    int print_points = atoi(argv[14]);
    int verbose_integration = atoi(argv[15]);
    size_t number_sample_points = static_cast<size_t>(atof(argv[16]));

    printf ("max_function_evals: %i\n"
            "min_function_evals: %i\n"
            "epsilon_abs: %g\n"
            "epsilon_rel: %g\n"
            "mode_tolerance: %g\n"
            "fraction_used_points: %g\n"
            "scale_increment: %g\n"
            "mass_fractions_file: %s\n"
            "posterior_output_file: %s\n"
            "points_output_file: %s\n"
            "pileup_input_file: %s\n"
            "selected_alg: %i\n"
            "collect_points: %i\n" 
            "print_points: %i\n"
            "verbose_integration: %i\n"
            "number_sample_points: %i\n",

            static_cast<int>(max_function_evals),
            static_cast<int>(min_function_evals),
            epsilon_abs,
            epsilon_rel,
            mode_tolerance,
            fraction_used_points,
            scale_increment,
            mass_fractions_file,
            posterior_output_file,
            points_output_file,
            pileup_input_file,
            selected_alg,
            collect_points,
            print_points,
            verbose_integration,
            static_cast<int>(number_sample_points));



    //lookup table to convert a pileup code to reduced code for
    //initial processing
    
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


    int last_index;
    int scanned;

    ErrorEstimate error_estimate;
    //error_estimate.Initialize(base_qual_prior_file);
   
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
        //int const* l = PileupSummary::base_to_index;

        //int sample_alphabet_size = strlen(sample_base_alphabet);

        //compile the posterior for this pileup summary
        //std::vector<SAMPLE> posterior_curve(sample_alphabet_size);

        std::vector<BaseComp> base_observations =
            ExpandTo4State(p._bases_upper, p._quality_codes, p._read_depth);

        BASE_COMP_COUNTS base_comp_counts = TallyBaseCompCounts(base_observations);


        //Find the mode point
        gsl_multimin_function minimizer_function;

        //Bounds full_bounds[] = { { 0.0, 1.0 }, { 0.0, 1.0 }, { 0.0, 1.0 } };

        minimizer_function.n = 3;
        minimizer_function.f = &gsl_simplex_hypercube_posterior;
        //minimizer_function.f = &gsl_simplex_posterior;

        
        bool add_teepee_component = true; //adds a very gently sloping teepee to the posterior to help find the peak

        const int initial_call_count = 0;

        PosteriorParams posterior_params = {
            &error_estimate, &base_comp_counts, 0.0, 
            static_cast<std::vector<WeightedSample> *>(NULL), 
            { { 0.0, 1.0 }, { 0.0, 1.0 }, { 0.0, 1.0 } },
            scale_increment,
            scale_increment,
            initial_call_count,
            0.0, //initial mode point
            { 0.0, 0.0, 0.0 } // mode point placeholder
        };
        

        minimizer_function.params = &posterior_params;
        gsl_vector *ss, *x, *x_cur, *x_prev, *x_jump;
        gsl_multimin_fminimizer * minimizer = 
            gsl_multimin_fminimizer_alloc(gsl_multimin_fminimizer_nmsimplex2, 3);

        x = gsl_vector_alloc(3);
        gsl_vector_set(x, 0, 0.25);
        gsl_vector_set(x, 1, 0.25);
        gsl_vector_set(x, 2, 0.25);

        ss = gsl_vector_alloc(3);
        gsl_vector_set_all(ss, 0.1);

        x_cur = gsl_vector_alloc(3);
        x_prev = gsl_vector_alloc(3);
        x_jump = gsl_vector_alloc(3);

        gsl_multimin_fminimizer_set(minimizer, &minimizer_function, x, ss);

       
        double prev_simplex_size = 10.0;
        double simplex_size = 1.0;
        double jump_size = 1.0;

        int max_multimin_iterations = 100000;
        int multimin_iteration = 0;
        while (gsl_multimin_test_size(simplex_size, mode_tolerance) == GSL_CONTINUE
               && multimin_iteration < max_multimin_iterations)
        {
            prev_simplex_size = simplex_size;
            gsl_multimin_fminimizer_iterate(minimizer);
            simplex_size = gsl_multimin_fminimizer_size(minimizer);
            ++multimin_iteration;
            if (multimin_iteration % 10000 == 0)
            {
//                 printf("simplex_size: %20.20g, iteration: %i\n",
//                        simplex_size, static_cast<int>(multimin_iteration));
            }
        }

//         printf("Found mode in %i iterations with simplex size %20.20f\n",
//                multimin_iteration, simplex_size);

        double mode[] = {
            gsl_vector_get(minimizer->x,0),
            gsl_vector_get(minimizer->x,1),
            gsl_vector_get(minimizer->x,2)
        };


        double expanded_mode[] = {
            gsl_vector_get(minimizer->x,0),
            gsl_vector_get(minimizer->x,1),
            gsl_vector_get(minimizer->x,2)
        };

        //error_estimate.expand_to_hypercube(mode, expanded_mode);
        error_estimate.contract_from_hypercube(expanded_mode, mode);

        //error_estimate.contract_from_hypercube(expanded_mode, mode);

        BaseComp mode_comp(mode[0], mode[1], mode[2], 1.0 - mode[0] - mode[1] - mode[2]);
        
        //deallocate minimizer and x vector
        gsl_multimin_fminimizer_free(minimizer);
        gsl_vector_free(x);
        gsl_vector_free(ss);

        posterior_params.best_log2_mode = error_estimate.Log2Posterior(base_comp_counts, mode_comp);

        //initially set the shift factor to the best log2 mode.
        posterior_params.shift_factor = posterior_params.best_log2_mode;

//         printf("Mode point: %f %f %f %f, Log2Posterior at mode = %f\n",
//                mode_comp.A, mode_comp.C, mode_comp.G, mode_comp.T, posterior_params.shift_factor);

        posterior_params.mode_point[0] = mode[0];
        posterior_params.mode_point[1] = mode[1];
        posterior_params.mode_point[2] = mode[2];

        
        int ndim = 3; //four nucleotides, sum to 1, implies three dimensions
        const int ncomp = 1; //only one component to the posterior function for now.

        int VERBOSE = verbose_integration ? 2 : 0;
        int LAST = 4;
        int MERSENNE = 8;

        std::vector<WeightedSample> weighted_samples[] = {
            std::vector<WeightedSample>(),
            std::vector<WeightedSample>(),
            std::vector<WeightedSample>(),
            std::vector<WeightedSample>()
        };


        char const* integrator_names[] = { "Vegas", "Suave", "Divonne", "Cuhre" };

        double integral[4][ncomp];
        double error[4][ncomp];
        double integral_correct_prob[4][ncomp];
        int nregions[4];
        int neval[4];
        int failcode[4];

        //bool RunMode[] = { true, true, true, true };
        bool RunMode[] = { false, false, false, false };
        RunMode[selected_alg] = true;

        int vegas_nstart = 100;
        int vegas_nincrease = 100;
        int suave_nnew = 10000;
        double suave_flatness = 1.0;

        int divonne_nregions, divonne_neval, divonne_failcode;
        
        int divonne_key1 = 11;
        int divonne_key2 = 11;
        int divonne_key3 = 1;
        
        int divonne_maxpass = 10;
        
        double divonne_border = 0.0;
        double divonne_maxchisq = 1e-30;
        double divonne_mindeviation = 1e-30;
        int divonne_ngiven = 1;
        int divonne_ldxgiven = ndim;
        double divonne_xgiven[3];
        int divonne_nextra = 0;

        int cuhre_key = 11; //a cubature rule!

        //use slice sampling to get a set of sample points
        int initial_sampling_range = 3 * 8;

         std::vector<BaseComp> sample_points = 
             SliceSampling(error_estimate, base_comp_counts, mode_comp, 
                           posterior_params.best_log2_mode,
                           initial_sampling_range,
                           number_sample_points);


        for (int algo = 0; algo != 4; ++algo)
        {

            posterior_params.samples = collect_points ? &weighted_samples[algo] : NULL;
            posterior_params.scale_factor = scale_increment;
            posterior_params.call_count = initial_call_count;

            if (RunMode[algo])
            {
                
                if (algo == 0)
                {
//                     Vegas(ndim, ncomp, cuba_posterior, &posterior_params, epsilon_rel, epsilon_abs, VERBOSE | LAST,
//                           min_function_evals, max_function_evals, 
//                           vegas_nstart, vegas_nincrease,
//                           &neval[algo], &failcode[algo], integral[algo], error[algo], integral_correct_prob[algo]);
                }
                else if (algo == 1)
                {
                    Suave(ndim, ncomp, cuba_posterior, &posterior_params, epsilon_rel, epsilon_abs, VERBOSE | LAST,
                          min_function_evals, max_function_evals, 
                          vegas_nstart, vegas_nincrease,
                          &nregions[algo], &neval[algo], &failcode[algo], 
                          integral[algo], error[algo], integral_correct_prob[algo]);
                }
                else if (algo == 2)
                {

                    REAL last_log2_mode = 0.0;

                    do
                    {

                        last_log2_mode = posterior_params.best_log2_mode;

                        posterior_params.scale_factor = 1.0;
                        posterior_params.shift_factor = posterior_params.best_log2_mode;
                        posterior_params.call_count = initial_call_count;

                        if (posterior_params.samples != NULL)
                        {
                            (*posterior_params.samples).clear();
                        }
                        error_estimate.expand_to_hypercube(posterior_params.mode_point, expanded_mode);
                        divonne_xgiven[0] = expanded_mode[0];
                        divonne_xgiven[1] = expanded_mode[1];
                        divonne_xgiven[2] = expanded_mode[2];

                        Divonne(ndim, ncomp, divonne_posterior, &posterior_params, 
                                epsilon_rel, epsilon_abs, LAST | MERSENNE | VERBOSE,
                                min_function_evals, max_function_evals, 
                                divonne_key1, divonne_key2, divonne_key3, divonne_maxpass,
                                divonne_border, divonne_maxchisq, divonne_mindeviation,
                                divonne_ngiven, divonne_ldxgiven, divonne_xgiven, divonne_nextra, NULL,
                                &nregions[algo], &neval[algo], &failcode[algo], 
                                integral[algo], error[algo], integral_correct_prob[algo]);
                        
                        printf("%i log2_mode: %20.20f, dif: %20.20f, (% .15f, % .15f, % .15f, % .15f), "
                               "%8s: %6.6g +- %6.6g (%6.6g)"
                               "%10i (reg) %10i (eval) %i (fail)\n",
                               summary._position,
                               posterior_params.best_log2_mode,
                               posterior_params.best_log2_mode - last_log2_mode,
                               posterior_params.mode_point[0],
                               posterior_params.mode_point[1],
                               posterior_params.mode_point[2],
                               1 
                               - posterior_params.mode_point[0]
                               - posterior_params.mode_point[1]
                               - posterior_params.mode_point[2],
                               integrator_names[algo], integral[algo][0], 
                               error[algo][0], integral_correct_prob[algo][0], 
                               nregions[algo], neval[algo], failcode[algo]);


                    }
                    while (last_log2_mode != posterior_params.best_log2_mode);

                }
                else if (algo == 3)
                {
//                     Cuhre(ndim, ncomp, cuhre_posterior, &posterior_params, epsilon_rel, epsilon_abs, VERBOSE | LAST,
//                           min_function_evals, max_function_evals, cuhre_key,
//                           &nregions[algo], &neval[algo], &failcode[algo], integral[algo], 
//                           error[algo], integral_correct_prob[algo]);
                }

//                 printf("%8s: %6.6g +- %6.6g (%6.6g)  %10i (reg) %10i (eval) %i (fail)\n",
//                        integrator_names[algo], integral[algo][0], error[algo][0], integral_correct_prob[algo][0], 
//                        nregions[algo], neval[algo], failcode[algo]);
            }
        }



        size_t print_limit = 0;

        for (size_t integrator = 0; integrator != 4; ++integrator)
        {
            if (print_limit > 0)
            {
                printf("\n**************************************************\n%s\n"
                       "%i samples\n"
                       "********************************************************\n",
                       integrator_names[integrator], 
                       static_cast<int>(weighted_samples[integrator].size()));
            }
            for (size_t s = 0; s != weighted_samples[integrator].size(); ++s)
            {
                if (s >= print_limit)
                {
                    break;
                }
                WeightedSample const& w = weighted_samples[integrator][s];
                printf("%s %7i %14.14f %14.14f %14.14f %14.14f: val: %14.14f wt: %14.14f\n",
                       integrator_names[integrator],
                       static_cast<int>(s), w.x[0], w.x[1], w.x[2], w.x[3], w.val, w.weight);
            }
        }

        //hack to select a subset of the samples
        for (size_t integrator = 0; integrator != 4; ++integrator)
        {
            std::vector<WeightedSample> & ws = weighted_samples[integrator];
            ws.erase(ws.begin(), ws.begin() + (1.0 - fraction_used_points) * ws.size());
            std::vector<WeightedSample>::iterator wsi;
            for (wsi = ws.begin(); wsi != ws.end(); ++wsi)
            {
                if ((*wsi).val == 0)
                {
                    //erase the zero-val element
                    //ws.erase(wsi);
                }
            }

        }

        
        //find desired marginal integral boundaries
        std::map<REAL, REAL> marginal_cdfs[4][4];

        for (size_t algo = 0; algo != 4; ++algo)
        {
            if (! RunMode[algo])
            {
                continue;
            }

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
                        "%s\t%g\t%i\t%7i\t%c", integrator_names[algo],
                        integral[algo][0], failcode[algo],
                        static_cast<int>(weighted_samples[algo].size()),
                        ErrorEstimate::nucleotides[base]);

                MarginalCumulativeDistribution(base, &weighted_samples[algo]);

                for (size_t mi = 0; mi != num_mass_fractions; ++mi)
                {

                    bool is_lower_bound = mass_fractions[mi] < 0.5;
                    REAL bound = 
                        FindIntegralBound(&weighted_samples[algo], base,
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

            fprintf(posterior_fh, 
                    "%s\t%g\t%i\t%7i\t%c", integrator_names[algo],
                    integral[algo][0], failcode[algo],
                    static_cast<int>(weighted_samples[algo].size()),
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
                        "%s\tlocus\tA\tC\tG\tT\tcA\tcC\tcG\tcT\t"
                        "integrand\tpoint_weight\n",
                        integrator_names[algo]);
                    
                for (size_t s = 0; s != weighted_samples[algo].size(); ++s)
                {
                        
                    WeightedSample & ws = weighted_samples[algo][s];
                    fprintf(points_fh, 
                            "%s\t%i\t%.10g\t%.10g\t%.10g\t%.10g"
                            "\t%.10g\t%.10g\t%.10g\t%.10g"
                            "\t%.10g\t%.10g\n",
                            integrator_names[algo],
                            summary._position, 
                            ws.x[0], ws.x[1], ws.x[2], ws.x[3], 
                            ws.cdf[0], ws.cdf[1], ws.cdf[2], ws.cdf[3],
                            ws.val, ws.weight);
                }
            }
        }
        
    }

    fclose(pileup_input_fh);
    fclose(posterior_fh);
    fclose(points_fh);

    return 0;
}
