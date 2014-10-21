#include <cstdio>
#include <cstdlib>

#include "pileup_tools.h"
#include "nucleotide_stats.h"
#include "samutil/file_utils.h"

int bqs_usage()
{
    fprintf(stderr,
            "\nTally {basecall, quality score, strand} counts in a sample pileup file\n"
            "Usage: dep bqs [options] sample.pileup sample.bqs\n"
            "Options:\n\n"
            "-t INT      number of threads to use [1]\n"
            "-m INT      number bytes of memory to use [100000000]\n"
            // "-F STRING   Fastq offset type if known (one of Sanger,Solexa,Illumina13,Illumina15) [None]\n"
            );
    return 1;
}

extern char *optarg;
extern int optind;


struct fastq_tally_input
{
    size_t thread_num;
    std::vector<char *>::iterator beg;
    std::vector<char *>::iterator end;
    size_t * counts; // do not own
    size_t min_quality_score;
    fastq_tally_input(size_t thread_num,
                      std::vector<char *>::iterator beg,
                      std::vector<char *>::iterator end,
                      size_t * counts,
                      size_t min_quality_score) :
        thread_num(thread_num), beg(beg), end(end), counts(counts),
        min_quality_score(min_quality_score)
    {
    }
    
};


void * fastq_tally_worker(void * args)
{
    fastq_tally_input * input = static_cast<fastq_tally_input *>(args);

    PileupSummary locus;
    std::vector<char *>::iterator it;
    for (it = input->beg; it != input->end; ++it)
    {
        locus.load_line(*it);
        locus.parse(input->min_quality_score);
        for (size_t c = 0; c != locus.counts.num_data; ++c)
        {
            input->counts[locus.counts.stats_index[c]] += locus.counts.raw_counts[c];
        }
    }
    pthread_exit((void*) 0);
}


int main_bqs(int argc, char ** argv)
{
    char c;
    size_t num_threads = 1;
    size_t max_mem = 100000000;
    // char const* fastq_type = "None";

    while ((c = getopt(argc, argv, "t:m:")) >= 0)
    {
        switch(c)
        {
        case 't': num_threads = static_cast<size_t>(atof(optarg)); break;
        case 'm': max_mem = static_cast<size_t>(atof(optarg)); break;
        // case 'F': fastq_type = optarg; break;
        default: return bqs_usage(); break;
        }
    }
    if (argc - optind != 2)
        return bqs_usage();

    char *pileup_input_file = argv[optind];
    char *bqs_output_file = argv[optind + 1];
    size_t min_quality_score = 0;

    //initialize fastq_type;

    size_t chunk_size = max_mem;
    char * chunk_buffer_in = new char[chunk_size + 1];

    int offset = fastq_offset(pileup_input_file, chunk_buffer_in, chunk_size);
    if (offset == -1)
    {
        fprintf(stderr, "dep bqs: Cannot continue.\n");
        return 1;
    }

    PileupSummary::set_offset(offset);
    
    FILE * pileup_input_fh = open_if_present(pileup_input_file, "r");
    FILE * bqs_output_fh = open_if_present(bqs_output_file, "w");
    
    
    char * last_fragment;
    size_t max_pileup_line_size = 1000000; // !!! fix this
    size_t bytes_read;
    size_t bytes_wanted = chunk_size - max_pileup_line_size;
    size_t * counts = new size_t[Nucleotide::num_bqs];
    size_t ** counts_t = new size_t *[num_threads];
    fastq_tally_input ** worker_input = new fastq_tally_input *[num_threads];
    pthread_t * threads = new pthread_t[num_threads];

    std::vector<char *> lines;

    for (size_t t = 0; t != num_threads; ++t)
    {
        counts_t[t] = new size_t[Nucleotide::num_bqs];
        std::fill(counts_t[t], counts_t[t] + Nucleotide::num_bqs, 0);
        worker_input[t] = new fastq_tally_input(t, lines.begin(), lines.end(), counts_t[t],
                                                min_quality_score);
    }

    size_t fread_nsec;
    size_t total_fread_nsec = 0;
    size_t total_bytes_read = 0;

    while (! feof(pileup_input_fh))
    {
        bytes_read = FileUtils::read_until_newline(chunk_buffer_in, bytes_wanted,
                                                   max_pileup_line_size, pileup_input_fh,
                                                   & fread_nsec);
        total_fread_nsec += fread_nsec;
        total_bytes_read += bytes_read;

        lines = FileUtils::find_complete_lines_nullify(chunk_buffer_in, & last_fragment);

        for (size_t t = 0; t != num_threads; ++t)
        {
            size_t b = range_chunk_offset(0, lines.size(), num_threads, t, true);
            size_t e = range_chunk_offset(0, lines.size(), num_threads, t, false);
            worker_input[t]->beg = lines.begin() + b;
            worker_input[t]->end = lines.begin() + e;

            int rc = pthread_create(&threads[t], NULL, &fastq_tally_worker, 
                                    static_cast<void *>(worker_input[t]));
            assert(rc == 0);
        }

        for (size_t t = 0; t != num_threads; ++t) {
            void *end;
            int rc = pthread_join(threads[t], &end);
            assert(0 == rc);
        }

    }
    fclose(pileup_input_fh);

    fprintf(stderr, "File reading metrics:  %Zu total bytes read in %Zu nanoseconds, %5.3f MB/s\n",
            total_bytes_read, total_fread_nsec,
            static_cast<float>(total_bytes_read) * 1000.0 / static_cast<float>(total_fread_nsec));

    delete chunk_buffer_in;
    delete threads;

    std::fill(counts, counts + Nucleotide::num_bqs, 0);
    for (size_t t = 0; t != num_threads; ++t)
    {
        for (size_t bqs = 0; bqs != Nucleotide::num_bqs; ++bqs)
        {
            counts[bqs] += worker_input[t]->counts[bqs];
        }
        delete counts_t[t];
        delete worker_input[t];
    }
    delete counts_t;
    delete worker_input;


    char basecall;
    size_t quality;
    size_t strand;
    for (size_t i = 0; i != Nucleotide::num_bqs; ++i)
    {
        Nucleotide::decode(i, &basecall, &quality, &strand);
        fprintf(bqs_output_fh,
                "%c\t%Zu\t%c\t%Zu\n",
                basecall,
                quality,
                (strand == Nucleotide::PLUS_STRAND) ? '+' : '-',
                counts[i]);
    }

    delete counts;

    close_if_present(bqs_output_fh);

    return 0;
}

