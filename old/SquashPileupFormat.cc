#include <cstdlib>
#include <cstdio>
#include <map>
#include <string>
#include <cstring>
#include <cmath>
#include <cassert>
#include <inttypes.h>

#include "henry/tools.h"
#include "henry/pileup_tools.h"

/*
  UNIX-style filter to convert samtools pileup format to squashed pileup format

  Squashed Pileup format
  "%i %i %c %i %i %i %i %i indel[ <indel counts from -n to n> ] 
  <position> <total_depth> <reference_base> <A_count> <C_count> <G_count> <T_count> N_count indel[ 

  34 48777 A 48201 258 106 124 86 indel[ 0 0 0 48775 1 2 2 ] seqs[ G CT TC GGGGG TCT]


  refposition refbase A C G T N i3 i2 i3 d1 d2 d3

*/

int main(int argc, char ** argv)
{

    int indel_histo_size = atoi(argv[1]);
    
    int last_index;
    int scanned;

    std::map<std::string, int>::iterator indel_iter;


    char sample_base_alphabet[10];
    
    char const* ansi_codes[] = { "", "", "", "" };
    //char const* ansi_codes[] = { "\033[0m", "\033[1m", "\033[22m", "\033[31m" };

    char const* ansi_close = ansi_codes[0];
    char const* ansi_encode = ansi_codes[1];
    char const* ansi_decode = ansi_codes[2];
    char const* ansi_red = ansi_codes[3];

    char char_values_fwd[] = "ACGTN";
    char char_values_rev[] = "acgtn";

    int scanned_fields;

    //fprintf(histogram_fh, "Position\tBhat\tChat\tD\tlogP\tlogG\n");

    while (! feof(stdin))
    {

        int next_char = fgetc(stdin);
        if (next_char == EOF)
        {
            break;
        }
        ungetc(next_char, stdin);

        PileupSummary summary(indel_histo_size);

        summary.load_line(stdin);
        
        PileupSummary const& p = summary;
        int const* l = PileupSummary::base_to_index;

        //make the majority base bold
        //make the majority base red if it mismatches the reference
        
        int max_count = 0, max_base = 0, current_count;
        for (int c = 0; char_values_fwd[c] != '\0'; ++c)
        {
            int nuc = static_cast<int>(char_values_fwd[c]);
            int rev = static_cast<int>(char_values_rev[c]);
            current_count = p.base_counts[static_cast<int>(l[nuc])] + 
                p.base_counts[static_cast<int>(l[rev])];

            max_count = std::max(current_count, max_count);
            max_base = current_count == max_count ? c : max_base;
        }

        bool is_variant = char_values_fwd[max_base] != toupper(summary._reference_base);
        char const* start_line = is_variant ? ansi_red : ansi_close;
        
        printf("%s", start_line);

        printf("%i %i %c ", p._position, p._read_depth, summary._reference_base);

        for (int c = 0; char_values_fwd[c] != '\0'; ++c)
        {
            bool is_reference_base = 
                char_values_fwd[c] == summary._reference_base ||
                char_values_rev[c] == summary._reference_base;

            int nuc = char_values_fwd[c];
            int rev = char_values_rev[c];

            char const* maybe_visual_open = is_reference_base ? ansi_encode : "";
            char const* maybe_close = 
                maybe_visual_open[0] != '\0' ? ansi_decode : "";

            printf("%s%i%s ",
                   maybe_visual_open,
                   p.base_counts[static_cast<int>(l[nuc])] + 
                   p.base_counts[static_cast<int>(l[rev])],
                   //p.base_qual_sums[l[nuc]] + p.base_qual_sums[l[rev]],
                   maybe_close);
        }
        printf("indel[ ");
        for (int i = -indel_histo_size; i <= indel_histo_size; ++i)
        {
            printf("%i ", summary.indel_counts[i]);
        }
        printf("] ");

        printf("seqs[");

//         int ansi_encode = 1;
//         int ansi_decode = 22;
        for (int i = -indel_histo_size; i <= indel_histo_size; ++i)
        {
            for (indel_iter = summary.indel_seqs[i].begin();
                 indel_iter != summary.indel_seqs[i].end(); ++indel_iter)
            {
                printf(" %s%s(%i)%s", ansi_encode, (*indel_iter).first.c_str(), 
                       (*indel_iter).second, ansi_decode);
            }
            if (i == 0) {
//                 ansi_encode = "\033[2m";
//                 ansi_decode = "\033[22m";
            }
        }

        printf("] ");
        printf("%s\n", ansi_close);

    }

    return 0;
}
