#include <cstdio>
#include <cstdlib>

#include "pileup_tools.h"
#include "error_estimate.h"
#include "samutil/file_utils.h"

struct base_qual_strand
{
    char basecall;
    size_t quality;
    char strand;
    bool operator<(base_qual_strand const& b) const
    {
        return this->basecall < b.basecall 
            || (this->basecall == b.basecall
                && (this->quality < b.quality
                    || (this->quality == b.quality
                        && (this->strand < b.strand))));
    }
                            
};


std::map<base_qual_strand, size_t>
load_base_qual_strand_data(char const* bases, 
                           char const* quality_codes,
                           size_t quality_code_offset,
                           size_t num_reads,
                           size_t min_quality_score)
{

    std::map<base_qual_strand, size_t> data;
    for (size_t r = 0; r != num_reads; ++r)
    {
        base_qual_strand bqs = 
            {
                static_cast<char>(toupper(bases[r])),
                QualityCodeToQuality(quality_codes[r], quality_code_offset),
                isupper(bases[r]) ? '+' : '-'
            };
        if (bqs.quality < min_quality_score)
        {
            continue;
        }

        if (data.find(bqs) == data.end())
        {
            data[bqs] = 0;
        }
        data[bqs]++;
    }
    return data;

}

int bqslocus_usage()
{
    fprintf(stderr, 
            "Usage: dep bqslocus input.pileup label base_qual_strand.cts base_strand.cts strand.cts\n\n"
            "Set any of base_qual_strand.cts, base_strand.cts, strand.cts to /dev/null\n"
            "if not needed.  (Will run faster as well.)\n\n"
            );
    return 1;
}


int main_bqslocus(int argc, char ** argv)
{

    if (argc != 5)
    {
        return bqslocus_usage();
    }

    char * pileup_input_file = argv[1];
    char * sample_string = argv[2];
    char * base_qual_strand_file = argv[3];
    char * base_strand_file = argv[4];
    char * strand_file = argv[5];

    FILE * pileup_input_fh = fopen(pileup_input_file, "r");
    if (pileup_input_fh == NULL)
    {
        fprintf(stderr, "Couldn't open pileup input file %s\n",
                pileup_input_file);
        exit(1);
    }

    FILE * base_qual_strand_fh;
    if (strcmp(base_qual_strand_file, "/dev/null") == 0)
    {
        base_qual_strand_fh = NULL;
    }
    else
    {
        base_qual_strand_fh = fopen(base_qual_strand_file, "w");
        if (base_qual_strand_fh == NULL)
        {
            fprintf(stderr, "Couldn't open base_qual_strand_file %s\n",
                    base_qual_strand_file);
        }
    }

    FILE * base_strand_fh;
    if (strcmp(base_strand_file, "/dev/null") == 0)
    {
        base_strand_fh = NULL;
    }
    else
    {
        base_strand_fh = fopen(base_strand_file, "w");
        if (base_strand_fh == NULL)
        {
            fprintf(stderr, "Couldn't open base_strand_file %s\n",
                    base_strand_file);
        }
    }


    FILE * strand_fh;
    if (strcmp(strand_file, "/dev/null") == 0)
    {
        strand_fh = NULL;
    }
    else
    {
        strand_fh = fopen(strand_file, "w");
        if (strand_fh == NULL)
        {
            fprintf(stderr, "Couldn't open strand_file %s\n",
                    strand_file);
        }
    }

    
    const int min_quality_score = 0;

    char line_label[1000];

    PileupSummary summary(0);

    size_t chunk_size = 1024 * 1024;
    char * chunk_buffer_in = new char[chunk_size + 1];

    FastqType ftype = summary.FastqFileType(pileup_input_file, chunk_buffer_in, chunk_size);

    PileupSummary::SetFtype(ftype);

    char const* nucleotides = "ACGT";

    size_t nbytes_read, nbytes_unused = 0;
    char * last_fragment;
    char * read_pointer = chunk_buffer_in;

    while (! feof(pileup_input_fh))
    {

        nbytes_read = fread(read_pointer, 1, chunk_size - nbytes_unused, pileup_input_fh);

        std::vector<char *> pileup_lines =
            FileUtils::find_complete_lines_nullify(chunk_buffer_in, & last_fragment);

        std::vector<char *>::iterator pit;
        read_pointer[nbytes_read] = '\0';

        //for (size_t l = 0; l != pileup_lines.size(); ++l)
        for (pit = pileup_lines.begin(); pit != pileup_lines.end(); ++pit)
        {

            bool succeeded = summary.load_line((*pit));

            if (! succeeded)
            {
                fprintf(stderr, "Couldn't parse pileup line\n");
                exit(1);
            }
        
            PileupSummary const& p = summary;
        
            std::map<base_qual_strand, size_t> base_qual_strand_data =
                load_base_qual_strand_data(p._bases_upper, p._quality_codes,
                                           PileupSummary::quality_code_offset,
                                           p._read_depth, min_quality_score);
        
        
            size_t effective_depth = base_qual_strand_data.size();

            sprintf(line_label, "%s\t%s\t%i\t%c\t%Zu\t%Zu", 
                    sample_string, summary._reference, 
                    summary._position, summary._reference_base, summary._read_depth,
                    effective_depth);

            std::map<std::pair<char, char>, size_t> base_strand_qualstats;
            std::map<std::pair<char, char>, size_t> base_strand_counts;

            std::map<char, size_t> strand_qualstats;
            std::map<char, size_t> strand_counts;

            //initialize base_strand_qualstats and base_strand_counts
            char const* base;
            char const* strand;
            char const* strands = "+-";

            for (strand = strands; *strand != '\0'; ++strand)
            {
                strand_counts[*strand] = 0;
                strand_qualstats[*strand] = 0;

                for (base = nucleotides; *base != '\0'; ++base)
                {
                    std::pair<char, char> strand_base(*base, *strand);
                    base_strand_qualstats[strand_base] = 0;
                    base_strand_counts[strand_base] = 0;
                }
            }


            for (std::map<base_qual_strand, size_t>::iterator bit = base_qual_strand_data.begin(); 
                 bit != base_qual_strand_data.end(); ++bit)
            {
                base_qual_strand const& datum = (*bit).first;
                size_t count = (*bit).second;

                char basecall = datum.basecall;
                char strand = datum.strand;

                std::pair<char, char> basecall_strand(basecall, strand);

                base_strand_qualstats[basecall_strand] += datum.quality * count;
                strand_qualstats[strand] += datum.quality * count;

                base_strand_counts[basecall_strand] += count;
                strand_counts[strand] += count;
            
                if (base_qual_strand_fh != NULL)
                {
                    fprintf(base_qual_strand_fh, "%s\t%c\t%c\t%Zu\t%Zu\n", line_label, 
                            basecall, strand, datum.quality, count);
                }
            }
        
            if (base_strand_fh != NULL)
            {
                std::map<std::pair<char, char>, size_t>::iterator cit;
                for (cit = base_strand_qualstats.begin(); cit != base_strand_qualstats.end(); ++cit)
                {
                    std::pair<char, char> basecall_strand = (*cit).first;
                    char basecall = basecall_strand.first;
                    char strand = basecall_strand.second;
                    size_t quality_sum = (*cit).second;

                    fprintf(base_strand_fh, "%s\t%c\t%c\t%Zu\t%Zu\n", line_label, 
                            basecall, strand, quality_sum, base_strand_counts[basecall_strand]);
                }
                fflush(base_strand_fh);
            }

            if (strand_fh != NULL)
            {
                std::map<char, size_t>::iterator cit;
                for (cit = strand_qualstats.begin(); cit != strand_qualstats.end(); ++cit)
                {
                    char strand = (*cit).first;
                    size_t quality_sum = (*cit).second;
                
                    fprintf(strand_fh, "%s\t%c\t%Zu\t%Zu\n", line_label, 
                            strand, quality_sum, strand_counts[strand]);
                }
                fflush(strand_fh);
            }

        }
        nbytes_unused = strlen(last_fragment);
        memmove(chunk_buffer_in, last_fragment, nbytes_unused);
        read_pointer = chunk_buffer_in + nbytes_unused;

    }
    
    if (base_qual_strand_fh != NULL)
    {
        fclose(base_qual_strand_fh);
    }
    if (base_strand_fh != NULL)
    {
        fclose(base_strand_fh);
    }
    if (strand_fh != NULL)
    {
        fclose(strand_fh);
    }

    return 0;

}
