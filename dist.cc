#include <fstream>
#include <algorithm>
#include <numeric>
#include <set>

#include <gsl/gsl_math.h>
#include <gsl/gsl_errno.h>

#include <string.h>

#include "run_comp_or_mode.h"
#include "comp_functor.h"
#include "tools.h"
#include "usage_strings.h"
#include "file_utils.h"
#include "dist_worker.h"
#include "sampling.h"
#include "pileup_tools.h"
#include "stats_tools.h"
#include "vcf.h"

#define MIN(a,b) ((a) < (b) ? (a) : (b))

int dist_usage()
{
    fprintf(stderr, 
            "\nUsage: dep dist [options] samples.rdb contig_order.rdb\n" 
            "Options:\n\n"
            "-d STRING   (optional) name of output distance file.  If absent, do not perform distance calculation.\n"
            "-c STRING   (optional) name of output composition file.  If absent, do not output composition marginal estimates.\n"
            "-V STRING   (optional) name of output vcf file.  If absent, do not output VCF file\n"
            "-i STRING   (optional) name of output indel distance file.  If absent, do not perform indel distance calculation.\n"
            "-S STRING   name of sample_pairings.rdb file (see description below).  Required if -d or -i are given.\n"
            "-y REAL     Minimum mutational distance to report as changed [0.2]\n"
            "-T INT      number of sample points used for tuning proposal distribution [1000]\n"
            "-x INT      number of sample points used for first-pass distance checking [100]\n"
            "-X REAL     quantile used to test first-pass distance [0.01]\n"
            "-f INT      number of sample points used for final quantiles estimation [1000]\n"
            "-q INT      minimum quality score to include bases as evidence [%s]\n"
            "-F STRING   Fastq offset type if known (one of Sanger, Solexa, Illumina13, Illumina15) [None]\n"
            "-g <empty>  if present, print extra pileup fields in the output distance file [absent]\n"

            "-t INT      number of threads to use [1]\n"
            "-m INT      number bytes of memory to use [%Zu]\n"
            "-v <empty>  if present, be verbose [absent]\n"
            "\n"
            "These options affect how sampling is done.\n"
            "\n"
            "-D STRING   distance quantiles file [\"0.005 0.05 0.5 0.95 0.995\\n\"]\n"
            "-C STRING   composition quantiles file [\"0.005 0.05 0.5 0.95 0.995\\n\"]\n"
            "-P INT      number of random sample point pairings to estimate Sampling/Sampling distance distribution [10000]\n"
            "-a INT      target autocorrelation offset.  once reached, proposal tuning halts [6]\n"
            "-I INT      maximum number of proposal tuning iterations [10]\n"
            "-M REAL     autocorrelation maximum value [0.2]\n"
            "-z REAL     gradient tolerance used in finding the posterior mode (gsl_multimin) [1e-5]\n"
            "-p FLOAT    alpha value for dirichlet prior [0.1]\n"
            "-l INT      maximum length of a pileup line in bytes.  {Nuisance parameter} [100000]\n"
            "\n"
            "samples.rdb has lines of <sample_id><tab></path/to/sample.jpd><tab></path/to/sample.pileup>\n"
            "sample_pairings.rdb has lines of <sample_id><tab><sample_id>\n"
            "defining which pairs of samples are to be compared.  <sample_id> correspond\n"
            "with those in samples.rdb\n"
            "\n"
            "contig_order.rdb has lines of <contig><tab><ordering>\n"
            "It must be consistent and complete with the orderings of the contigs mentioned in\n"
            "all pileup input files\n"
            ,
            Usage::quality_string,
            1024 * 1024 * 1024 * 4l // 4 GB memory
            );
    return 1;
}

extern char *optarg;
extern int optind;



int main_dist(int argc, char **argv)
{

    size_t num_threads = 1;
    size_t max_mem = 1024l * 1024l * 1024l * 4l;

    struct posterior_settings pset = {
        1e-5,   // gradient_tolerance
        3000,   // max_modefinding_iterations
        10,     // max_tuning_iterations
        1000,   // tuning_num_points
        1000,   // final_num_points
        6,      // autocor_max_offset
        0.2,    // autocor_max_value
        30,     // initial_autocor_offset
        6,      // target_autocor_offset
        true,   // is_log_integrand
        62 * 3  // initial_sampling_range
    };


    size_t prelim_num_points = 1e2;
    float prelim_quantile = 0.01;
    size_t num_sample_point_pairs = 1e4;

    size_t max_pileup_line_size = 1e6;

    int print_pileup_fields = 0;

    float prior_alpha = 0.1;

    size_t min_quality_score = 5;

    const char *dist_quantiles_file = NULL;
    const char *comp_quantiles_file = NULL;

    const char *dist_file = NULL;
    const char *comp_file = NULL;
    const char *vcf_file = NULL;
    const char *indel_dist_file = NULL;
    const char *fastq_type = NULL;

    const char *sample_pairings_file = NULL;

    float min_dist_to_report = 0.2;

    bool verbose = false;

    char c;
    while ((c = getopt(argc, argv, "d:c:V:i:S:y:t:m:T:x:X:f:gP:a:I:M:z:p:q:F:vD:C:l:")) >= 0)
    {
        switch(c)
        {
        case 'd': dist_file = optarg; break;
        case 'c': comp_file = optarg; break;
        case 'V': vcf_file = optarg; break;
        case 'i': indel_dist_file = optarg; break;
        case 'S': sample_pairings_file = optarg; break;
        case 'y': min_dist_to_report = atof(optarg); break;
        case 't': num_threads = static_cast<size_t>(atof(optarg)); break;
        case 'm': max_mem = static_cast<size_t>(atof(optarg)); break;
        case 'T': pset.tuning_num_points = static_cast<size_t>(atof(optarg)); break;
        case 'x': prelim_num_points = static_cast<size_t>(atof(optarg)); break;
        case 'X': prelim_quantile = atof(optarg); break;
        case 'f': pset.final_num_points = static_cast<size_t>(atof(optarg)); break;
        case 'P': num_sample_point_pairs = static_cast<size_t>(atof(optarg)); break;
        case 'a': pset.target_autocor_offset = static_cast<size_t>(atof(optarg)); break;
        case 'I': pset.max_tuning_iterations = static_cast<size_t>(atof(optarg)); break;
        case 'M': pset.autocor_max_value = atof(optarg); break;
        case 'z': pset.gradient_tolerance = atof(optarg); break;
	case 'p': prior_alpha = atof(optarg); break;
        case 'q': min_quality_score = static_cast<size_t>(atoi(optarg)); break;
        case 'F': fastq_type = optarg; break;
        case 'v': verbose = true; break;
        case 'D': dist_quantiles_file = optarg; break;
        case 'C': comp_quantiles_file = optarg; break;
        case 'g': print_pileup_fields = 1; break;
        case 'l': max_pileup_line_size = static_cast<size_t>(atof(optarg)); break;
        default: return dist_usage(); break;
        }
    }
    if (argc - optind != 2)
    {
        return dist_usage();
    }

    const char *samples_file = argv[optind];
    const char *contig_order_file = argv[optind + 1];

    // consistency checks
    if (! sample_pairings_file && (dist_file || indel_dist_file))
    {
        // either provide dist_file and sample_pairings_file, or none at all
        fprintf(stderr, "If -d or -i are given, you must also provide -S.\n");
        exit(5);
    }

    if (! dist_file && ! comp_file && ! indel_dist_file)
    {
        fprintf(stderr, "Error: You must provide at least one of -d or -c or -i.  "
                "Otherwise, there is nothing to calculate\n");
        exit(5);
    }

    gsl_set_error_handler_off();

    std::vector<char *> sample_label;
    std::vector<char *> jpd_input_file;
    std::vector<char *> pileup_input_file;

    FILE *samples_fh = open_if_present(samples_file, "r");
    while (! feof(samples_fh))
    {
        char *jpd = new char[1000];
        char *pileup = new char[1000];
        char *label = new char[100];
        fscanf(samples_fh, "%s\t%s\t%s\n", label, jpd, pileup);
        sample_label.push_back(label);
        jpd_input_file.push_back(jpd);
        pileup_input_file.push_back(pileup);
    }
    fclose(samples_fh);

    size_t num_samples = pileup_input_file.size();
    size_t chunk_size = 1024 * 1024 * 64;
    char *chunk_buffer_in = new char[chunk_size + 1];

    int offset;

    if (fastq_type)
        offset = fastq_type_to_offset(fastq_type);
    else
        offset = fastq_offset(pileup_input_file[0], chunk_buffer_in, chunk_size);

    if (offset == -1)
    {
        fprintf(stderr, "Could not determine fastq type of this pileup file.\n");
        return 1;
    }

    PileupSummary::set_offset(offset);

    delete chunk_buffer_in;

    FILE **pileup_input_fh = new FILE *[num_samples];
    for (size_t s = 0; s != num_samples; ++s)
    {
        pileup_input_fh[s] = open_if_present(pileup_input_file[s], "r");
    }

    FILE *dist_fh = open_if_present(dist_file, "w");
    FILE *comp_fh = open_if_present(comp_file, "w");
    FILE *vcf_fh = open_if_present(vcf_file, "w");
    FILE *indel_dist_fh = open_if_present(indel_dist_file, "w");

    double *dist_quantiles;
    double *comp_quantiles;
    double default_quantiles[] = { 0.005, 0.05, 0.5, 0.95, 0.995 };

    size_t num_dist_quantiles;
    size_t num_comp_quantiles;

    if (dist_quantiles_file == NULL)
    {
        dist_quantiles = new double[5];
        std::copy(default_quantiles, default_quantiles + 5, dist_quantiles);
        num_dist_quantiles = 5;
    }
    else
    {
        dist_quantiles = ParseNumbersFile(dist_quantiles_file, & num_dist_quantiles);
    }
    if (comp_quantiles_file == NULL)
    {
        comp_quantiles = new double[5];
        std::copy(default_quantiles, default_quantiles + 5, comp_quantiles);
        num_comp_quantiles = 5;
    }
    else
    {
        comp_quantiles = ParseNumbersFile(comp_quantiles_file, & num_comp_quantiles);
    }

    double prior_alphas[NUM_NUCS];
    for (size_t i = 0; i != NUM_NUCS; ++i)
        prior_alphas[i] = prior_alpha;

    // parse pairings file
    size_t num_pairings, *pair_sample1, *pair_sample2;

    {
        FILE *sample_pairings_fh = open_if_present(sample_pairings_file, "r");
        std::vector<size_t> pair1, pair2;
        std::vector<size_t> *trg[] = { & pair1, & pair2 };
        char label[2][100];
        while (! feof(sample_pairings_fh))
        {
            std::vector<char *>::iterator it;
            fscanf(sample_pairings_fh, "%s\t%s\n", label[0], label[1]);
            
            for (size_t p = 0; p != 2; ++p)
            {
                for (it = sample_label.begin(); it != sample_label.end(); ++it)
                {
                    if (strcmp(*it, label[p]) == 0)
                    {
                        (*trg[p]).push_back(std::distance(sample_label.begin(), it));
                        break;
                    }
                }
                if (it == sample_label.end())
                {
                    fprintf(stderr, "Error: label %s in %s not found in samples file %s\n",
                            label[p], sample_pairings_file, samples_file);
                    exit(10);
                }
            }
        }
        fclose(sample_pairings_fh);
        
        num_pairings = pair1.size();
        pair_sample1 = new size_t[num_pairings];
        pair_sample2 = new size_t[num_pairings];
        std::copy(pair1.begin(), pair1.end(), pair_sample1);
        std::copy(pair2.begin(), pair2.end(), pair_sample2);
    }
    
    std::map<char const*, size_t, ltstr> contig_order;
    // parse contig_order file
    {
        FILE *contig_order_fh = open_if_present(contig_order_file, "r");
        while (! feof(contig_order_fh))
        {
            char *contig = new char[100];
            size_t order;
            fscanf(contig_order_fh, "%s\t%zu\n", contig, & order);
            contig_order.insert(std::make_pair(contig, order));
        }
        fclose(contig_order_fh);
    }


    std::vector<std::vector<char *> > pileup_lines(num_samples);
    std::vector<std::vector<char *>::iterator> current(num_samples);
    std::vector<std::vector<char *>::iterator> bound(num_samples);
    char **chunk_buffer = new char *[num_samples];
    size_t *chunk_read_size = new size_t[num_samples];

    size_t single_buffer_size = max_mem / (num_samples * 2); // approximate...
    size_t max_dist_line_size = distance_quantiles_locus_bytes(num_dist_quantiles);
    size_t max_comp_line_size = marginal_quantiles_locus_bytes(num_comp_quantiles);
    size_t max_vcf_line_size = vcf_locus_bytes(num_samples);

    if (single_buffer_size < max_pileup_line_size)
    {
        fprintf(stderr, "Error: single_buffer_size = %Zu, which is smaller than max_pileup_line_size %Zu.  "
                "Either decrease option -l or increase option -m\n",
                single_buffer_size, max_pileup_line_size);
        exit(10);
    }


    // initialize
    for (size_t s = 0; s != num_samples; ++s)
    {
        current[s] = pileup_lines[s].begin();
        bound[s] = pileup_lines[s].end();
        chunk_buffer[s] = new char[single_buffer_size + 1];
        chunk_read_size[s] = 0;
    }
    std::vector<char *>::iterator least_upper_bound = bound[0]; // set to arbitrary 'least'

    dist_worker_input **worker_inputs = new dist_worker_input*[num_threads];
    
    // initialize all generic structures in the worker_inputs
    for (size_t t = 0; t != num_threads; ++t)
    {
        worker_inputs[t] = 
            new dist_worker_input(t, num_samples, num_pairings, num_sample_point_pairs,
                                  dist_quantiles, num_dist_quantiles,
                                  comp_quantiles, num_comp_quantiles,
                                  min_dist_to_report,
                                  prelim_num_points,
                                  prelim_quantile,
                                  pset.final_num_points,
                                  print_pileup_fields,
                                  NULL,
                                  NULL,
                                  NULL,
                                  NULL,
                                  pair_sample1,
                                  pair_sample2,
                                  & contig_order);
        
        for (size_t s = 0; s != num_samples; ++s)
        {
            worker_inputs[t]->worker[s] = 
                new posterior_wrapper(jpd_input_file[s],
                                      prior_alphas,
                                      min_quality_score,
                                      comp_quantiles,
                                      num_comp_quantiles,
                                      sample_label[s],
                                      NULL,
                                      NULL,
                                      pset,
                                      verbose);
        }
    }

    char *last_fragment;
    less_locus_position less_locus;
    equal_locus_position equal_locus;
    less_locus.contig_order = &contig_order;
    equal_locus.contig_order = &contig_order;

    /* The main chunk loop.  At the beginning of each iteration, the
      range [current[s], bound[s]) represents the set of input lines
      in sample S that have been read from file but not yet processed.
      If the range is empty, more of the file is read and the range is
      updated.
      
     */
    size_t move_threshold = single_buffer_size / 10;
    size_t bytes_wanted = single_buffer_size - max_pileup_line_size - 1;
    bool all_input_read = false;
    size_t total_pass_bytes_read;
    
    size_t fread_nsec; // number of nanoseconds to read this chunk.
    size_t total_fread_nsec = 0;
    size_t total_bytes_read = 0;

    // print the vcf header
    if (vcf_fh != NULL)
    {
        fprintf(vcf_fh,
               "##fileformat=VCFv4.1\n"
               "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT");

        for (size_t s = 0; s != num_samples; ++s)
        {
            fprintf(vcf_fh, "\t%s", sample_label[s]);
        }
        fprintf(vcf_fh, "\n");
    }
    fflush(vcf_fh);

    while (! all_input_read)
    {
        // for each sample, determine whether the supply of lines is
        // exhausted, if so, refresh it.
        total_pass_bytes_read = 0;

        for (size_t s = 0; s != num_samples; ++s)
        {
            bool do_full_reload = current[s] == pileup_lines[s].end();
            size_t bytes_read = 0;
            if (do_full_reload)
            {
                // do a full reload
                bytes_read = FileUtils::read_until_newline(chunk_buffer[s], bytes_wanted,
                                                           max_pileup_line_size, pileup_input_fh[s],
                                                           & fread_nsec);
                
                total_fread_nsec += fread_nsec;
                chunk_read_size[s] = bytes_read;

                assert(bytes_read == 0 || chunk_buffer[s][chunk_read_size[s] - 1] == '\n');

                pileup_lines[s] = 
                    FileUtils::find_complete_lines_nullify(chunk_buffer[s], &last_fragment);

                assert(chunk_buffer[s][chunk_read_size[s] - 1] == '\0');

                current[s] = pileup_lines[s].begin();
                bound[s] = pileup_lines[s].end();
            }
            else
            {
                size_t left_shift = *current[s] - chunk_buffer[s];
                if (left_shift > move_threshold)
                {
                    // do a shift.
                    assert(chunk_buffer[s][chunk_read_size[s] - 1] == '\0');
                    size_t new_read_offset = chunk_read_size[s] - left_shift;

                    memmove(chunk_buffer[s], chunk_buffer[s] + left_shift, new_read_offset);

                    // update all moved pointers
                    for (std::vector<char *>::iterator lit = current[s]; lit != pileup_lines[s].end(); ++lit)
                    {
                        *lit -= left_shift;
                    }
                    bytes_read = FileUtils::read_until_newline(chunk_buffer[s] + new_read_offset,
                                                               bytes_wanted - new_read_offset, 
                                                               max_pileup_line_size, pileup_input_fh[s],
                                                               & fread_nsec);

                    total_fread_nsec += fread_nsec;

                    chunk_read_size[s] = new_read_offset + bytes_read;

                    assert(bytes_read == 0 || chunk_buffer[s][chunk_read_size[s] - 1] == '\n');

                    assert(chunk_read_size[s] < single_buffer_size);

                    // need to splice the vector
                    std::vector<char *> tmp =
                        FileUtils::find_complete_lines_nullify(chunk_buffer[s] + new_read_offset, &last_fragment);

                    assert(chunk_buffer[s][chunk_read_size[s] - 1] == '\0');

                    pileup_lines[s].erase(pileup_lines[s].begin(), current[s]);
                    pileup_lines[s].insert(pileup_lines[s].end(), tmp.begin(), tmp.end());
                    current[s] = pileup_lines[s].begin();
                    bound[s] = pileup_lines[s].end();
                }
                else
                {
                    // do nothing
                }
            }
            total_pass_bytes_read += bytes_read;
        }
        all_input_read = (total_pass_bytes_read == 0);

        total_bytes_read += total_pass_bytes_read;

        size_t max_num_loci = 0; // number of loci in one sample
        size_t sample_with_max = 0;

        // if any of pileup_lines[s] are non-empty, global_s will be
        // set to the one whose rbegin() has the least genomic position
        size_t global_s = 0;
        
        for (size_t s = 0; s != num_samples; ++s)
        {
            if ((! pileup_lines[s].empty()) 
                && (pileup_lines[global_s].empty()
                    || less_locus(*pileup_lines[s].rbegin(), *pileup_lines[global_s].rbegin())
                    ))
            {
                global_s = s;
            }
        }
        
        // refresh all upper bounds
        // this ensures that no pairs are processed until it is known
        // that the required loci are parsed
        for (size_t s = 0; s != num_samples; ++s)
        {
            if (all_input_read || pileup_lines[global_s].empty())
            {
                bound[s] = pileup_lines[s].end();
            }
            else
            {
                bound[s] = 
                    std::upper_bound(current[s], pileup_lines[s].end(), 
                                     *pileup_lines[global_s].rbegin(), less_locus);
            }

            size_t num_loci = std::distance(current[s], bound[s]);
            max_num_loci = std::max(max_num_loci, num_loci);
            sample_with_max = num_loci == max_num_loci ? s : sample_with_max;
        }

        // prepare workers.  divide up each interval according to the
        // interval with the most loci.
        pthread_t *threads = new pthread_t[num_threads];

        for (size_t t = 0; t != num_threads; ++t)
        {
            // find the guide positions for the maximal sample
            // guide_end is guaranteed to be a valid iterator except when t == num_threads - 1
            // in that case, it equals bound[sample_with_max]
            std::vector<char *>::iterator guide_end = 
                current[sample_with_max] + (((t + 1) * max_num_loci) / num_threads);

            size_t max_sample_loci = 0;
            size_t total_loci = 0;
            for (size_t s = 0; s != num_samples; ++s)
            {
                worker_inputs[t]->beg[s] = (t == 0) ? current[s] : worker_inputs[t-1]->end[s];
                worker_inputs[t]->end[s] = (t == num_threads - 1) 
                    ? bound[s]
                    : std::upper_bound(current[s], bound[s], *guide_end, less_locus);

                size_t num_loci = std::distance(worker_inputs[t]->beg[s], worker_inputs[t]->end[s]);
                total_loci += num_loci;
                max_sample_loci = std::max(num_loci, max_sample_loci);
            }

            // pessimistically, assume that there is zero overlap between any loci.
            // therefore, generate one dist line for every locus position and pairing
            size_t dist_output_size = max_dist_line_size * (total_loci * num_pairings);
            size_t comp_output_size = max_comp_line_size * (max_sample_loci * num_samples);
            size_t vcf_output_size = max_vcf_line_size * (max_sample_loci * num_samples);
            size_t indel_dist_output_size = dist_output_size;

            worker_inputs[t]->out_dist = (dist_fh != NULL) ? new char[dist_output_size + 1] : NULL;
            worker_inputs[t]->out_comp = (comp_fh != NULL) ? new char[comp_output_size + 1] : NULL;
            worker_inputs[t]->out_vcf = (vcf_fh != NULL) ? new char[vcf_output_size + 1] : NULL;
            worker_inputs[t]->out_indel_dist = (indel_dist_fh != NULL) ? new char[indel_dist_output_size + 1] : NULL;

            // dispatch threads
            int rc = pthread_create(&threads[t], NULL, &dist_worker, static_cast<void *>(worker_inputs[t]));
            assert(rc == 0);
        }


        for (size_t t = 0; t < num_threads; ++t) {
            int rc = pthread_join(threads[t], NULL);
            assert(0 == rc);
        }

        // write output

        if (dist_fh != NULL)
        {
            for (size_t t = 0; t < num_threads; ++t) {
                fwrite(worker_inputs[t]->out_dist, 1, strlen(worker_inputs[t]->out_dist), dist_fh);
                delete[] worker_inputs[t]->out_dist;
            }        
            fflush(dist_fh);
        }
        if (comp_fh != NULL)
        {
            for (size_t t = 0; t < num_threads; ++t) {
                fwrite(worker_inputs[t]->out_comp, 1, strlen(worker_inputs[t]->out_comp), comp_fh);
                delete[] worker_inputs[t]->out_comp;
            }        
            fflush(comp_fh);
        }
        if (vcf_fh != NULL)
        {
            for (size_t t = 0; t < num_threads; ++t) {
                fwrite(worker_inputs[t]->out_vcf, 1, strlen(worker_inputs[t]->out_vcf), vcf_fh);
                delete[] worker_inputs[t]->out_vcf;
            }        
            fflush(vcf_fh);
        }
        if (indel_dist_fh != NULL)
        {
            for (size_t t = 0; t < num_threads; ++t) {
                fwrite(worker_inputs[t]->out_indel_dist, 1, strlen(worker_inputs[t]->out_indel_dist), indel_dist_fh);
                delete[] worker_inputs[t]->out_indel_dist;
            }        
            fflush(indel_dist_fh);
        }

        delete[] threads;

        // update the ranges to reflect the set of loci processed
        for (size_t s = 0; s != num_samples; ++s)
        {
            current[s] = bound[s];
        }
    }

    fprintf(stderr, "File reading metrics:  %Zu total bytes read in %Zu nanoseconds, %5.3f MB/s\n",
            total_bytes_read, total_fread_nsec,
            static_cast<float>(total_bytes_read) * 1000.0 / static_cast<float>(total_fread_nsec));

    if (dist_fh != NULL) { fclose(dist_fh); }
    if (comp_fh != NULL) { fclose(comp_fh); }
    if (vcf_fh != NULL) { fclose(vcf_fh); }

    // cleanup
    for (size_t s = 0; s != num_samples; ++s)
    {
        fclose(pileup_input_fh[s]);
        delete[] chunk_buffer[s];
        delete[] jpd_input_file[s];
        delete[] pileup_input_file[s];
        delete[] sample_label[s];
    }
    for (size_t t = 0; t != num_threads; ++t)
    {
        // this owns its workers, so no need to delete them individually
        delete worker_inputs[t];
    }

    for (std::map<char const*, size_t, ltstr>::iterator cit = contig_order.begin(); 
         cit != contig_order.end(); ++cit)
    {
        delete (*cit).first;
    }
    delete[] pair_sample1;
    delete[] pair_sample2;

    delete[] chunk_buffer;
    delete[] pileup_input_fh;
    delete[] dist_quantiles;
    delete[] comp_quantiles;
    delete[] chunk_read_size;

    return 0;
}
