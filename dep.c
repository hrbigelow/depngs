#include "dep.h"

#include <stdio.h>
#include <string.h>

int usage()
{
    fprintf(stderr, "\n"
            "dep (diversity estimation probabilistically)\n"
            "Author: Henry Bigelow (hrbigelow@gmail.com)\n\n"
            "Usage:\n\n"
            "dep dist       Sample base composition and pairwise base differences with confidence intervals\n"
            "dep pileup     Multi-way pileup\n"
            "dep pstats     Output quality-binned basecall distribution locus histogram\n"
            // "dep diststats  Compute a table of pairwise 4D Dirichlet distance classifications\n"
            // "dep mergestats Merge 2 or more diststats files\n"
            // "dep simp       Simulate pileup file\n"
            // "dep simc       Simulate loci compositions\n"
            // "dep bqslocus   Tally {basecall, quality score, strand} counts per locus\n"
            // "dep bqs        Tally {basecall, quality score, strand} counts overall\n"
            // "dep bqs2jpd    Expand {basecall, quality score, strand} counts to jpd\n"
            // "dep pug        Pileup Grep.  grep a list of loci from a large pileup file using binary search\n"
            "\n"
            );
    return 1;
}

int main(int argc, char **argv)
{
    if (argc < 2) 
        return usage();

    else if (strcmp(argv[1], "dist") == 0)
        return main_dist(argc - 1, argv + 1);

    else if (strcmp(argv[1], "pileup") == 0)
        return main_pileup(argc - 1, argv + 1);

    else if (strcmp(argv[1], "pstats") == 0)
        return main_pstats(argc - 1, argv + 1);
    
    // else if (strcmp(argv[1], "diststats") == 0)
    //     return main_diststats(argc - 1, argv + 1);

    // else if (strcmp(argv[1], "mergestats") == 0)
    //     return main_mergestats(argc - 1, argv + 1);

    // else if (strcmp(argv[1], "simp") == 0)
    //     return main_simp(argc - 1, argv + 1);

    // else if (strcmp(argv[1], "simc") == 0)
    //     return main_simc(argc - 1, argv + 1);

    // else if (strcmp(argv[1], "bqs") == 0)
    //     return main_bqs(argc - 1, argv + 1);

    // else if (strcmp(argv[1], "bqslocus") == 0)
    //     return main_bqslocus(argc - 1, argv + 1);

    // else if (strcmp(argv[1], "bqs2jpd") == 0)
    //     return main_bqs2jpd(argc - 1, argv + 1);

    /* else if (strcmp(argv[1], "pug") == 0) */
    /*     return main_pug(argc - 1, argv + 1); */

    else
    {
        fprintf(stderr, "Error: unrecognized command '%s'\n", argv[1]);
        return 2;
    }
    return 0;
}
