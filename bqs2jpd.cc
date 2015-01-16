#include "tools.h"

extern "C" {
#include "nucleotide_stats.h"
}

int bqs2jpd_usage()
{
    fprintf(stderr,
            "\nCompute joint probability distribution P(basecall, quality score, strand) from a set of .bqs counts\n"
            "\nUsage: dep bqs2jpd sample.bqs sample.jpd\n\n");
    return 1;
}

int main_bqs2jpd(int argc, char ** argv)
{

    if (argc != 3)
    {
        return bqs2jpd_usage();
    }

    char basecall;
    size_t quality;
    char strand;
    float count;
    float error_prob;
    float data_prob[4];
    size_t basecall_index;


    FILE * counts_fh = fopen(argv[1], "r");

    char const* jpd_output_file = argv[2];
    FILE * jpd_output_fh = open_if_present(jpd_output_file, "w");

    while (! feof(counts_fh))
    {
        fscanf(counts_fh, "%c\t%zi\t%c\t%f\n", &basecall, &quality, &strand, &count);
        error_prob = QualityToErrorProb(quality);

        std::fill(data_prob, data_prob + 4, (error_prob / 3.0) * count);
        basecall_index = base_to_index(basecall);
        data_prob[basecall_index] = (1.0 - error_prob) * count;
        
        fprintf(jpd_output_fh, "%c_%Zu_%c\t%8.6g\t%8.6g\t%8.6g\t%8.6g\n",
                basecall, quality, strand, data_prob[0], data_prob[1],
                data_prob[2], data_prob[3]);
    }

    fclose(counts_fh);
    close_if_present(jpd_output_fh);

    return 0;
}
