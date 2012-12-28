#include <cstdlib>
#include <cstdio>
#include <ctime>

#include <numeric>
#include <algorithm>

#include <gsl/gsl_rng.h>
#include <gsl/gsl_math.h>


#include <sys/timeb.h>

#include "henry/sampling.h"
#include "henry/error_estimate.h"
#include "henry/tools.h"
#include "henry/simulation.h"

/*

Generate, line by line, a simulated pileup file.
Each line is produced as a generated set of N bases

Desired output:

pileup-formatted, generated samples, one per line
corresponding generative parameter file consisting of:
mgen  Mgen  msam  Msam  A C G T

*/




int main(int argc, char ** argv)
{

    char * real_composition_file = argv[1];
    size_t sample_size = static_cast<size_t>(atof(argv[2]));
    char * base_qual_prior_file = argv[3];
//     char * parameters_file = argv[4];

    double real_base_composition[10000];

    size_t num_compositions = 
        ParseNumbersFile(real_composition_file, real_base_composition)
        / 4;

    char * founder_bases = new char[sample_size];
    char * called_bases = new char[sample_size];
    char * called_bases_printed = new char[sample_size];
    int * measured_quality = new int[sample_size];
    
//     FILE * parameters_fh = fopen(parameters_file, "w");
//     if (parameters_fh == NULL)
//     {
//         fprintf(stderr, "Couldn't open parameters output file %s for writing\n",
//                 parameters_file);
//         exit(1);
//     }

    gsl_rng * rand_gen = gsl_rng_alloc(gsl_rng_taus);
    timeb millitime;
    ftime(& millitime);
    gsl_rng_set(rand_gen, millitime.millitm);

    int pileup_line = 0;

    char reference[200];

    BASE_COMP_COUNTS empty_counts;

    ErrorEstimate error_estimate;

    BaseQualStrandReader reader;
    NucleotideStats prior_data = reader.read_from_rdb(base_qual_prior_file);
    error_estimate.set_prior_data(&prior_data);

    size_t const nbases = 4;
    size_t const nquals = ErrorEstimate::NQUALS;
    size_t const nstrands = 2;

    size_t const nraw_data = nbases * nquals * nstrands;

    

    //initialize all data probabilities P(I|b).  For each b, there are 34 x 4 P(I|b).
    double data_probability[nbases * nraw_data];
    double data_probability_sum;

    for (size_t ci = 0; ci != num_compositions; ++ci)
    {
        double * this_composition = real_base_composition + (ci * 4);

        for (size_t founder_bi = 0; founder_bi != nbases; ++founder_bi)
        {
            for (size_t major_bi = 0; major_bi != nbases; ++major_bi)
            {
                for (size_t qual = 0; qual != nquals; ++qual)
                {
                    for (size_t strand_num = 0; strand_num != 2; ++strand_num)
                    {
                        size_t flat_index = 
                            (founder_bi * nraw_data)
                            + (major_bi * nquals * nstrands)
                            + (strand_num * nquals)
                            + qual;

                        DNAStrand strand = strand_num == 0 ? POS_STRAND : NEG_STRAND;
                               
                        BaseComp image_comp(major_bi, qual, strand);
                        
                        data_probability[flat_index] = 
                            error_estimate.data_probability(image_comp, founder_bi)
                            * this_composition[founder_bi];
                    }
                }
            }
            
        }

        //this should sum to 1
        data_probability_sum =
            std::accumulate(data_probability, 
                            data_probability + (nbases * nbases * nquals), 0.0);

        assert(gsl_fcmp(data_probability_sum, 1.0, 1e-10) == 0);
        
        size_t majority_comp_index = 
            std::distance(this_composition, 
                          std::max_element(this_composition, 
                                           this_composition + 4));

        char majority_comp = ErrorEstimate::nucleotides[majority_comp_index];


        for (size_t si = 0; si != sample_size; ++si)
        {

            size_t flat_index =
                SampleDiscreteDistribution(rand_gen, data_probability,
                                           data_probability_sum);

            size_t founder_base_index = flat_index / nraw_data;
            size_t major_base_index = (flat_index % nraw_data) / (nquals * nstrands);
            size_t strand = (flat_index % nraw_data) / nquals;
            size_t quality_score = flat_index % (nquals * nstrands);

            founder_bases[si] = ErrorEstimate::nucleotides[founder_base_index];
            char called_base = ErrorEstimate::nucleotides[major_base_index];

            called_bases[si] = called_base;
            called_bases_printed[si] = 
                majority_comp == called_base ? '.' : called_base;
                
            measured_quality[si] = quality_score;
        }

        //print out a line of pileup

        float founder_base_comp[nbases];
        for (size_t bi = 0; bi != nbases; ++bi)
        {
            founder_base_comp[bi] =
                static_cast<float>(std::count(founder_bases, 
                                              founder_bases + sample_size,
                                              ErrorEstimate::nucleotides[bi])) 
                / static_cast<float>(sample_size);
        }

        sprintf(reference, "%05.5f_%05.5f_%05.5f_%05.5f_"
                "%05.5f_%05.5f_%05.5f_%05.5f",
                this_composition[0],
                this_composition[1],
                this_composition[2],
                this_composition[3],
                founder_base_comp[0],
                founder_base_comp[1],
                founder_base_comp[2],
                founder_base_comp[3]);

        PrintPileupLine(stdout, reference, pileup_line++,
                        majority_comp, called_bases_printed, 
                        measured_quality, sample_size);

    }
    
//     fclose(parameters_fh);
    delete founder_bases;
    delete called_bases;
    delete called_bases_printed;
    delete measured_quality;
    gsl_rng_free(rand_gen);

}
