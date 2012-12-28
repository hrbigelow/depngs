#include <cstdio>
#include <cstdlib>

#include "henry/pileup_tools.h"
#include "henry/tools.h"

//UNIX-style filter to accept a pileup format, and output a very
//expanded RDB representation of each base.


//Pileup format Input:  reference  position  reference_base  read_depth  pileup_codes  qualities
//Output format      :  reference  position  reference_base  sample_base  sample_quality



int main()
{
    PileupSummary::initialize();
    PileupSummary summary(0);
    while (! feof(stdin))
    {
        if (! summary.load_line(stdin))
        {
            exit(1);
        }
        for (int i = 0; i != summary._read_depth; ++i)
        {
            PileupSummary const& s = summary;
            printf("%s\t%i\t%c\t%c\t%i\n", s._reference,
                   s._position, s._reference_base,
                   s._bases[i], QualityCodeToQuality(s._quality_codes[i]));
        }
        fscanf(stdin, " "); //eat white space.  will this result in EOF being reached?

    }
    return 0;
}

       
