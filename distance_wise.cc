#include <cstdlib>
#include <cstddef>
#include <cstdint>
#include <cstdio>
#include <algorithm>
#include <gsl/gsl_math.h>
#include <pthread.h>
#include <cassert>
#include <vector>
#include <string.h>
#include <samutil/file_utils.h>

/*
  load an RDB file, and efficiently perform the pairwise query:

  select 
  a.snp_seq_changes_id,
  a.snp_sample_id sample1_id,
  b.snp_sample_id sample2_id, 
  sqrt(sum(power(a.comp - b.comp, 2))) dist_to_modes,
  sqrt(sum(power(case when a.
  from snp_mode_indiv a, snp_mode_indiv b, snp_pairwise_comp c
  where a.snp_seq_changes_id = b.snp_seq_changes_id
  and a.snp_sample_id = c.snp_sample1_id
  and b.snp_sample_id = c.snp_sample2_id
  group by a.snp_seq_changes_id, c.comparison_id
 */


/*
  input should be organized by (top to bottom), with cardinality given ()

  snp_sample_id (35)
  locus (11,000,000)
  hyp_base (4)
  stats (4) (rank, mode, lo, hi)

  output should be:

  comparison_id (N)
  locus (11,000,000)
  stats (3)

 */

struct input_stat
{
    uint32_t rank; // to ensure good alignment
    float mode;
    float lo;
    float hi;
};


struct output_stat
{
    uint16_t sample1_id;
    uint16_t sample2_id;
    uint32_t locus_id;
    float dist_mode;
    float dist_lo;
    float dist_hi;
};

struct calc_input
{
    size_t locus_id;
    size_t sample1_id;
    size_t sample2_id;
    input_stat * sample1_beg;
    input_stat * sample1_end;
    input_stat * sample2_beg;
    output_stat * out_start;
};


struct comparison
{
    size_t sample1_id;
    size_t sample2_id;
};

inline float ssd(float a1, float a2, float a3, float a4,
                 float b1, float b2, float b3, float b4)
{
    return 
        gsl_pow_2(a1 - b1)
        + gsl_pow_2(a2 - b2)
        + gsl_pow_2(a3 - b3)
        + gsl_pow_2(a4 - b4);
}


void * calc_distances(void * args)
{
    calc_input * input = static_cast<calc_input *>(args);
    input_stat * sample1 = input->sample1_beg;
    input_stat * sample2 = input->sample2_beg;
    output_stat * out_ptr = input->out_start;

    size_t num_bases = 4;
    size_t locus_id = input->locus_id;
    while (sample1 != input->sample1_end)
    {
            
        (*out_ptr).sample1_id = input->sample1_id;
        (*out_ptr).sample2_id = input->sample2_id;
        (*out_ptr).locus_id = locus_id;
        if ((*sample1).rank == 5 || (*sample2).rank == 5)
        {
            // one or both of the input combinations are unknown.
            // set metrics to dummy metrics
            (*out_ptr).dist_mode = -1.0;
            (*out_ptr).dist_lo = 0.0;
            (*out_ptr).dist_hi = 1000.0;
        }
        else
        {
            (*out_ptr).dist_mode =
                sqrt(
                     ssd(
                         (*sample1).mode, (*(sample1 + 1)).mode, (*(sample1 + 2)).mode, (*(sample1 + 3)).mode,
                         (*sample2).mode, (*(sample2 + 1)).mode, (*(sample2 + 2)).mode, (*(sample2 + 3)).mode)
                     );
                
            // this is an estimate of the minimum distance these two loci could be apart
            // this is based on the notion that, if the principle base (rank 0) is 
            float s_up[2][4];
            float s_dn[2][4];
            input_stat * stat;
            for (size_t s = 0; s != 2; ++s)
            {
                for (size_t b = 0; b != 4; ++b)
                {
                    stat = ((s == 0) ? sample1 : sample2) + b;
                    if ((*stat).rank == 0)
                    {
                        s_up[s][b] = (*stat).lo;
                        s_dn[s][b] = (*stat).hi;
                    }
                    else if ((*stat).rank == 1)
                    {
                        s_up[s][b] = (*stat).hi;
                        s_dn[s][b] = (*stat).lo;
                    }
                    else
                    {
                        s_up[s][b] = (*stat).mode;
                        s_dn[s][b] = (*stat).mode;
                    }
                }
            }
                
            // try all four combinations of distances
            float d[4];
            d[0] = ssd(s_up[0][0], s_up[0][1], s_up[0][2], s_up[0][3],
                       s_up[1][0], s_up[1][1], s_up[1][2], s_up[1][3]);
                
            d[1] = ssd(s_up[0][0], s_up[0][1], s_up[0][2], s_up[0][3],
                       s_dn[1][0], s_dn[1][1], s_dn[1][2], s_dn[1][3]);
                
            d[2] = ssd(s_dn[0][0], s_dn[0][1], s_dn[0][2], s_dn[0][3],
                       s_up[1][0], s_up[1][1], s_up[1][2], s_up[1][3]);
                
            d[3] = ssd(s_dn[0][0], s_dn[0][1], s_dn[0][2], s_dn[0][3],
                       s_dn[1][0], s_dn[1][1], s_dn[1][2], s_dn[1][3]);
                
            // find the max and min distances
            (*out_ptr).dist_lo = 100000.0;
            (*out_ptr).dist_hi = 0.0;
            for (size_t di = 0; di != 4; ++di)
            {
                (*out_ptr).dist_lo = std::min((*out_ptr).dist_lo, d[di]);
                (*out_ptr).dist_hi = std::max((*out_ptr).dist_hi, d[di]);
            }
            (*out_ptr).dist_lo = sqrt((*out_ptr).dist_lo);
            (*out_ptr).dist_hi = sqrt((*out_ptr).dist_hi);
        }   
        sample1 += num_bases;
        sample2 += num_bases;
        ++locus_id;
        ++out_ptr;
    }
    pthread_exit((void *) 0);
}


struct line_range
{
    size_t sample_input_size;
    std::vector<char *>::iterator start;
    std::vector<char *>::iterator end;
    input_stat * out;
};


void * init_input(void * args)
{
    line_range * input = static_cast<line_range *>(args);
    std::vector<char *>::iterator cur_line = input->start;
    size_t sample_input_size = input->sample_input_size;

    input_stat * stat;
    size_t sample_id;
    size_t locus_id;
    char hyp_base;
    size_t rank;
    float mode, lo, hi;
    size_t num_bases = 4;

    while (cur_line != input->end)
    {
        sscanf(*cur_line, "%zi\t%zi\t%c\t%zi\t%f\t%f\t%f",
               & locus_id,
               & sample_id,
               & hyp_base,
               & rank,
               & mode,
               & lo,
               & hi);

        stat = input->out + (sample_input_size * sample_id + (num_bases * locus_id));
        switch (hyp_base)
        {
        case 'A': break;
        case 'C': stat += 1; break;
        case 'G': stat += 2; break;
        case 'T': stat += 3; break;
        }
        (*stat).rank = rank;
        (*stat).mode = mode;
        (*stat).lo = lo;
        (*stat).hi = hi;

        ++cur_line;
    }
    pthread_exit((void *) 0);
}



int main(int argc, char ** argv)
{
    size_t num_threads = static_cast<size_t>(atof(argv[1]));
    size_t num_samples = static_cast<size_t>(atof(argv[2]));
    size_t num_loci = static_cast<size_t>(atof(argv[3]));
    size_t chunk_size = static_cast<size_t>(atof(argv[4]));
    char * data_file = argv[5];
    char * comp_file = argv[6];

    size_t num_bases = 4;

    size_t total_input_size = num_samples * num_loci * num_bases;
    size_t sample_input_size = total_input_size / num_samples;

    input_stat * input_stats = new input_stat[total_input_size];
    input_stat * in_ptr = input_stats;
    input_stat * in_end = input_stats + total_input_size;

    size_t num_comparisons = 0;
    comparison * comparisons = new comparison[1000];
    comparison * comp_ptr = comparisons;
    // parse comparisons file
    FILE * comp_fh = fopen(comp_file, "r");
    if (comp_fh == NULL)
    {
        fprintf(stderr, "Error, couldn't open file %s\n", comp_file);
        exit(1);
    }

    while (! feof(comp_fh))
    {
        fscanf(comp_fh, "%zu\t%zu\n", & (*comp_ptr).sample1_id, & (*comp_ptr).sample2_id);
        comp_ptr++;
        num_comparisons++;
    }
    fclose(comp_fh);


    output_stat * output_stats = new output_stat[num_loci];
    output_stat * out_ptr = output_stats;
    output_stat * out_end = output_stats + num_loci;
    
    // since we are doing 'perfect insertion sort', we don't know which values will be missing.
    // therefore, initialize everything to the default of 'missing'
    while (in_ptr != in_end)
    {
        (*in_ptr++).rank = 5;
    }

    // load input
    FILE * in_fh = fopen(data_file, "r");
    if (in_fh == NULL)
    {
        fprintf(stderr, "Error, couldn't open data file %s\n", data_file);
        exit(1);
    }
    size_t max_line_length = 100;
    size_t nbytes_read;
    size_t nbytes_want = chunk_size - max_line_length;
    char * chunk_buffer = new char[chunk_size];
    char * dummy; // needed for 'find_complete_lines_nullify'

    while (! feof(in_fh))
    {
        nbytes_read = fread(chunk_buffer, 1, nbytes_want, in_fh);
        chunk_buffer[nbytes_read] = '\0';

        if (nbytes_read == nbytes_want)
        {
            // keep reading until we hit the next newline
            char * cur = chunk_buffer + nbytes_read;
            char * test = fgets(cur, max_line_length, in_fh);
            if (test != cur)
            {
                fprintf(stderr, "Error: unexpected line size\n");
                exit(2);
            }
            nbytes_read += strlen(cur);
        }
        else
        {
            // didn't get the requested bytes.  should be at the end of the file...
            if (chunk_buffer[nbytes_read - 1] != '\n')
            {
                fprintf(stderr, "Error: File doesn't end in a newline\n");
                exit(1);
            }
            else
            {
                // all good.  do nothing
            }
        }

        std::vector<char *> lines =
            FileUtils::find_complete_lines_nullify(chunk_buffer, & dummy);

        line_range * line_ranges = new line_range[num_threads];
        size_t work_load = lines.size() / num_threads;
        for (size_t t = 0; t != num_threads; ++t)
        {
            line_ranges[t].sample_input_size = sample_input_size;
            line_ranges[t].start = lines.begin() + (t * work_load);
            line_ranges[t].end = (t == (num_threads - 1)) ? lines.end() : lines.begin() + ((t + 1) * work_load);
            line_ranges[t].out = input_stats;
        }

        pthread_t * threads = new pthread_t[num_threads];
        for (size_t t = 0; t != num_threads; ++t)
        {
            int rc = pthread_create(&threads[t], NULL,
                                    & init_input,
                                    static_cast<void *>(& line_ranges[t]));
            assert(rc == 0);
        }

        for (size_t t = 0; t < num_threads; ++t) {
            int rc = pthread_join(threads[t], NULL);
            assert(0 == rc);
        }

        delete[] threads;
        delete[] line_ranges;

    }

    delete chunk_buffer;
    fclose(in_fh);

    // compute output stats
    input_stat * sample1;
    input_stat * sample2;
    input_stat * sample1_end;

    size_t max_out_line_size = 50;
    char * chunk_buffer_out = new char[max_out_line_size * num_loci];

    for (size_t c = 0; c != num_comparisons; ++c)
    {
        size_t sample1_id = comparisons[c].sample1_id;
        size_t sample2_id = comparisons[c].sample2_id;
        sample1 = input_stats + (sample_input_size * sample1_id);
        sample2 = input_stats + (sample_input_size * sample2_id);
        sample1_end = sample1 + sample_input_size;
        out_ptr = output_stats;

        // here, instead of a single loop, num_threads loops
        calc_input * calc_inputs = new calc_input[num_threads];
        size_t work_load = sample_input_size / (num_threads * num_bases);
        for (size_t t = 0; t != num_threads; ++t)
        {
            calc_inputs[t].sample1_id = sample1_id;
            calc_inputs[t].sample2_id = sample2_id;
            calc_inputs[t].locus_id = t * work_load;
            calc_inputs[t].sample1_beg = sample1 + (t * work_load * num_bases);
            calc_inputs[t].sample1_end = (t == (num_threads - 1)) ? sample1_end : sample1 + ((t + 1) * work_load * num_bases);
            calc_inputs[t].sample2_beg = sample2 + (t * work_load * num_bases);
            calc_inputs[t].out_start = out_ptr + (t * work_load);
        }

        pthread_t * threads = new pthread_t[num_threads];
        for (size_t t = 0; t != num_threads; ++t)
        {
            int rc = pthread_create(&threads[t], NULL,
                                    & calc_distances,
                                    static_cast<void *>(& calc_inputs[t]));
            assert(rc == 0);
        }

        for (size_t t = 0; t < num_threads; ++t) {
            int rc = pthread_join(threads[t], NULL);
            assert(rc == 0);
        }

        delete[] threads;
        delete[] calc_inputs;

        // write output stats
        out_ptr = output_stats;
        char * out_buf_ptr = chunk_buffer_out;
        while (out_ptr != out_end)
        {
            out_buf_ptr += 
            sprintf(out_buf_ptr,
                    "%u\t%u\t%u\t%6.5f\t%6.5f\t%6.5f\n",
                    (*out_ptr).sample1_id,
                    (*out_ptr).sample2_id,
                    (*out_ptr).locus_id,
                    (*out_ptr).dist_mode,
                    (*out_ptr).dist_lo,
                    (*out_ptr).dist_hi);
            ++out_ptr;
        }
        fwrite(chunk_buffer_out, 1, (out_buf_ptr - chunk_buffer_out), stdout);
        fflush(stdout);
        
    }
    
    delete chunk_buffer_out;
    delete input_stats;
    delete output_stats;
    delete comparisons;
}
