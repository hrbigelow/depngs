#include <gsl/gsl_rng.h>
#include <gsl/gsl_math.h>
#include <sys/timeb.h>

int simp_usage()
{
    fprintf(stderr, 
            "Usage: dep simp FTYPE input.jpd dist_string max_depth n_rep sim1.pileup sim2.pileup dist.rdb\n"
            "FTYPE is one of 'Sanger (+33)' or 'Solexa (+64)'\n\n"
            "Simulate a pair of pileup files, for each combination of (depth1, depth2, distance, rep)\n"
            "D1 and D2 vary in [0, max_depth], distance is chosen from dist_string, and\n"
            "rep varies in [0, n_rep)\n"
            "\n"
            "The dist.rdb file has fields of '<position><tab><distance>'\n"
            );
    return 1;
}


struct pair_comp {
    float c1[4];
    float c2[4];
};


struct base_qual {
    char call, qual;
    unsigned strand;
};

/* Generate a base from a normalized composition */
unsigned gen_base(float *comp, gsl_rng *rg)
{
    unsigned i = 0;
    float q = gsl_rng_uniform(rg), cc = comp[0];
    while (q > cc) cc += comp[++i];
    return i;
}

#define MAX_DIST 100
unsigned qoffset;
float dist[MAX_DIST];
unsigned n_dist;

void init_dist(const char *csv)
{
    n_dist = 0;
    float *d = dist;
    char *loc;
    while (csv != '\0' && n_dist != MAX_DIST)
    {
        if (*csv == ',') ++csv;
        *d = strtod(csv, &loc);
        csv = loc;
        ++d;
        ++n_dist;
    }
}


static struct nucleotide_stats stats;

struct base_qual gen_call(unsigned bi, gsl_rng *rg)
{
    size_t code, quality, strand;

    unsigned i = 0;
    float q = gsl_rng_uniform(rg), code = cpd[bi][0];
    while (q > code) code += cpd[bi][++i];

    struct base_qual bq;
    decode_nucleotide(code, &bq.call, &quality, &strand);
    bq.qual = (char)(quality + qoffset);
    bq.strand = (unsigned)strand;
}

int main_simp(int argc, char **argv)
{
    if (argc != 9) return simp_usage();

    char *ftype_string = argv[1];
    char *base_qual_params_file = argv[2];
    char *dist_string = argv[3];
    unsigned max_depth = strtod(argv[4]);
    unsigned n_rep = strtod(argv[5]);
    char *sim1_file = argv[6];
    char *sim2_file = argv[7];
    char *dist_file = argv[8];

    if (strcmp(ftype_string, "Sanger") == 0) qoffset = 33;
    else if (strcmp(ftype_string, "Solexa") == 0) qoffset = 64;
    else fprintf(stderr, "Error: FTYPE must be 'Sanger' or 'Solexa'\n");

    FILE *sim1_fh = fopen(sim1_file, "w");
    if (! sim1_fh)
    {
        fprintf(stderr, "Couldn't open sim1_file %s for writing\n", sim1_file);
        exit(1);
    }

    FILE *sim2_fh = fopen(sim2_file, "w");
    if (! sim2_fh)
    {
        fprintf(stderr, "Couldn't open sim2_file %s for writing\n", sim2_file);
        exit(1);
    }

    FILE *dist_fh = fopen(dist_file, "w");
    if (! dist_fh)
    {
        fprintf(stderr, "Couldn't open dist_file %s for writing\n", dist_file);
        exit(1);
    }


    nucleotide_stats_initialize(base_qual_params_file, &stats);


    init_dist(dist_string);

    gsl_rng *rg = gsl_rng_alloc(gsl_rng_taus);
    timeb millitime;
    ftime(& millitime);
    gsl_rng_set(rg, millitime.millitm);

    model_params.initialize(base_qual_params_file);

    unsigned b, p = 1;
    struct pair_comp pc;
    struct base_qual bq;

    static char nucs[] = "ACGT";

    unsigned max_d1, max_d2, d1, d2, r, di;
    for (d1_max = 0; d1_max != max_depth; ++d1_max)
        for (d2_max = 0; d2_max != max_depth; ++d2_max)
            for (r = 0; r != n_rep; ++r)
                for (di = 0; di != n_dist; ++di)
                {
                    /* generate two random points from the simplex at
                       the specified distance. */
                    pc = gen_pair_comp(dist[di]);
                    
                    /* simulate from both, to the specified depth */
                    for (s = 0; s != 2; ++s)
                    {
                        float *c = s == 0 ? &pc.c1 : &pc.c2;
                        unsigned d_max = s == 0 ? d1_max : d2_max;
                        FILE **fhp = s == 0 ? sam1_fh : sam2_fh;

                        for (d = 0; d != d_max; ++d)
                        {
                            b = gen_base(c, rg);
                            bq = gen_call(b, rg);
                            calls[d] = b.strand == NUC_PLUS_STRAND 
                                ? b.call
                                : tolower(b.call);

                            quals[d] = b.qual;
                        }
                        calls[d_max] = '\0';
                        quals[d_max] = '\0';

                        fprintf(*fhp, "chr1\t%u\t%c\t%s\t%s\n",
                                p, 'N', d_max, calls, quals);
                    }
                    fprintf(dist_fh, "%u\t%g\n", p, dist[di]);
                    ++p;
                }

    gsl_rng_free(rg);

    fclose(sim1_fh);
    fclose(sim2_fh);
    fclose(dist_fh);

    return 0;
}
