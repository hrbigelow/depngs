#include <cstdio>
#include <cstdlib>

#include "pileup_tools.h"
#include "base_qual_strand_reader.h"
#include "nucleotide_stats.h"
#include "samutil/file_utils.h"

int bqs_usage()
{
    fprintf(stderr,
            "\nTally {basecall, quality score, strand} counts in a sample pileup file\n"
            "Usage: dep bqs sample.pileup sample.bqs\n\n");
    return 1;
}

int main_bqs(int argc, char ** argv)
{

    if (argc != 3)
    {
        return bqs_usage();
    }

    char * pileup_input_file = argv[1];
    char * bqs_output_file = argv[2];
    size_t const MAX_QUALITY = 255;
    size_t min_quality_score = 0;

    char const* strands = "+-";

    //initialize fastq_type;
    PileupSummary pileup(0);

    size_t chunk_size = 1024 * 1024;
    char * chunk_buffer_in = new char[chunk_size + 1];

    FastqType ftype = pileup.FastqFileType(pileup_input_file, chunk_buffer_in, chunk_size);

    if (ftype == None)
    {
        fprintf(stderr, "Error: Couldn't determine quality scale for pileup input file %s\n",
                pileup_input_file);
        exit(1);
    }

    PileupSummary::SetFtype(ftype);
    
    BaseQualStrandReader reader;
    // reader.initialize(pileup_input_file);

    FILE * pileup_input_fh = open_if_present(pileup_input_file, "r");
    FILE * bqs_output_fh = open_if_present(bqs_output_file, "w");
    
    double fake_nuc_frequency[] = { 1, 1, 1, 1 };

    //all we need is a name mapping.
    size_t datum_index = 0;

    JPD_DATA fake_jpd;
    for (size_t b = 0; b != 4; ++b)
    {
        for (size_t q = 0; q <= MAX_QUALITY; ++q)
        {
            for (size_t s = 0; s != 2; ++s)
            {
                BaseQualStrandReader::Datum datum = 
                    { Nucleotide::bases_upper[b], q, strands[s], datum_index, 0 };

                fake_jpd.insert(std::make_pair(datum.name(),nuc_frequency(fake_nuc_frequency)));
                ++datum_index;
            }
        }
    }

    NucleotideStats fake_stats(fake_jpd.size());
    fake_stats.initialize(fake_jpd);

    double * all_counts = new double[fake_jpd.size()];
    std::fill(all_counts, all_counts + fake_jpd.size(), 0.0);

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
            LocusSummary locus = 
                reader.get_next_locus(fake_stats, (*pit),
                                      static_cast<void const*>(& min_quality_score));
            
            for (size_t raw_index = 0; raw_index != locus.num_distinct_data; ++raw_index)
            {
                size_t data_index = locus.stats_index[raw_index];
                all_counts[data_index] += locus.raw_counts[raw_index];
            }
        }
        nbytes_unused = strlen(last_fragment);
        memmove(chunk_buffer_in, last_fragment, nbytes_unused);
        read_pointer = chunk_buffer_in + nbytes_unused;
    }
    fclose(pileup_input_fh);

    delete chunk_buffer_in;

    size_t max_quality_present = 0;

    std::map<std::string, size_t>::const_iterator datum_iter;

    //find the maximum nonzero quality present
    for (datum_iter = fake_stats.name_mapping.begin();
         datum_iter != fake_stats.name_mapping.end();
         ++datum_iter)
    {
        BaseQualStrandReader::Datum datum = 
            reader.get_datum_from_name((*datum_iter).first);
        size_t datum_index = (*datum_iter).second;

        if (all_counts[datum_index] > 0.0)
        {
            max_quality_present = std::max(max_quality_present, datum.quality);
        }
        
    }

    for (datum_iter = fake_stats.name_mapping.begin();
         datum_iter != fake_stats.name_mapping.end();
         ++datum_iter)
    {
        BaseQualStrandReader::Datum datum = 
            reader.get_datum_from_name((*datum_iter).first);
        size_t datum_index = (*datum_iter).second;

        if (datum.quality <= max_quality_present)
        {
            fprintf(bqs_output_fh, "%c\t%Zu\t%c\t%Zu\n",
                    datum.called_base,
                    datum.quality,
                    datum.strand,
                    static_cast<size_t>(all_counts[datum_index]));
        }
    }
    
    delete all_counts;
    close_if_present(bqs_output_fh);

    return 0;
}

