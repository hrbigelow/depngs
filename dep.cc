#include "dep.h"

#include <cstdio>
#include <cstring>

int usage()
{
    fprintf(stderr, "\n"
            "dep (diversity estimation probabilistically)\n"
            "Author: Henry Bigelow (hrbigelow@gmail.com)\n\n"
            "Usage:\n\n"
            "dep comp       Estimate per-locus base composition of diverse population\n"
            "dep mode       Estimate most probable per-locus base composition\n"
            "dep dist       Estimate sample-pairwise base composition distance with error bars\n"
            // "dep discomp    Estimate per-locus discrete genotype of clonal population\n"
            // "dep anomaly    Scores the degree of data anomaly\n"
            "dep simp       Simulate pileup file\n"
            "dep simc       Simulate loci compositions\n"
            "dep bqslocus   Tally {basecall, quality score, strand} counts per locus\n"
            "dep bqs        Tally {basecall, quality score, strand} counts overall\n"
            "dep bqs2jpd    Expand {basecall, quality score, strand} counts to jpd\n"
            "dep pug        Pileup Grep.  grep a list of loci from a large pileup file using binary search.\n"
            "\n"
            );
    return 1;
}

int main(int argc, char *argv[])
{

    if (argc < 2)
    {
        return usage();
    }
    else if (strcmp(argv[1], "comp") == 0)
    {
        return main_comp(argc - 1, argv + 1);
    }
    else if (strcmp(argv[1], "mode") == 0)
    {
        return main_mode(argc - 1, argv + 1);
    }
    else if (strcmp(argv[1], "dist") == 0)
    {
        return main_dist(argc - 1, argv + 1);
    }
    // else if (strcmp(argv[1], "discomp") == 0)
    // {
    //     return main_discomp(argc - 1, argv + 1);
    // }
    // else if (strcmp(argv[1], "anomaly") == 0)
    // {
    //     return main_anomaly(argc - 1, argv + 1);
    // }
    else if (strcmp(argv[1], "simp") == 0)
    {
        return main_simp(argc - 1, argv + 1);
    }
    else if (strcmp(argv[1], "simc") == 0)
    {
        return main_simc(argc - 1, argv + 1);
    }
    else if (strcmp(argv[1], "bqs") == 0)
    {
        return main_bqs(argc - 1, argv + 1);
    }
    else if (strcmp(argv[1], "bqslocus") == 0)
    {
        return main_bqslocus(argc - 1, argv + 1);
    }
    else if (strcmp(argv[1], "bqs2jpd") == 0)
    {
        return main_bqs2jpd(argc - 1, argv + 1);
    }
    else if (strcmp(argv[1], "pug") == 0)
    {
        return main_pug(argc - 1, argv + 1);
    }
    else
    {
        fprintf(stderr, "Error: unrecognized command '%s'\n", argv[1]);
        return 2;
    }
    return 0;
}
