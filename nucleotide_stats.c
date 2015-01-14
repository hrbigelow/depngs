#include "nucleotide_stats.h"

#include <stdio.h>
#include <string.h>
#include <assert.h>

struct nucleotide_stats make_nucleotide_stats()
{
    struct nucleotide_stats s;
    size_t b;
    for (b = 0; b != 4; ++b)
    {
        s.complete_jpd[b] = s.jpd_buffer + (b * NUC_NUM_BQS);
        s.founder_base_likelihood[b] = s.cpd_buffer + (b * NUC_NUM_BQS);
    }
    return s;
}


void nucleotide_stats_initialize(const char *rdb_file, struct nucleotide_stats *s)
{
    FILE *rdb_fh = fopen(rdb_file, "r");
    if (rdb_fh == NULL)
    {
        fprintf(stderr, "Couldn't open data jpd file %s\n", rdb_file);
        exit(5);
    }
    memset(s->jpd_buffer, 0, sizeof(s->jpd_buffer));
    double counts[4];

    char basecall, strand;
    int quality;
    size_t index_code, b, bi, di, fi;

    while (! feof(rdb_fh))
    {
        (void)fscanf(rdb_fh, "%c_%i_%c\t%lf\t%lf\t%lf\t%lf\n", 
                     &basecall, &quality, &strand, 
                     counts, counts+1, counts+2, counts+3);
        
        for (bi = 0; bi != 4; ++bi)
            if (counts[bi] < 0)
            {
                fprintf(stderr, "%s: found negative count for %c_%i_%c.\n", 
                        __func__, basecall, quality, strand);
                exit(11);
            }
        if (counts[0] + counts[1] + counts[2] + counts[3] == 0) continue;
        index_code = 
            encode_nucleotide(basecall, quality, 
                              (strand == '+' ? NUC_PLUS_STRAND : NUC_MINUS_STRAND));

        for (b = 0; b != 4; ++b)
            s->complete_jpd[b][index_code] = counts[b];
    }
    fclose(rdb_fh);

    double jpd_sum = 0, jpd_inv_sum;
    for (fi = 0; fi != NUC_NUM_FBQS; ++fi) jpd_sum += s->jpd_buffer[fi];
    jpd_inv_sum = 1.0 / jpd_sum;
    for (fi = 0; fi != NUC_NUM_FBQS; ++fi) s->jpd_buffer[fi] *= jpd_inv_sum;

    for (b = 0; b != 4; ++b)
    {
        s->founder_base_marginal[b] = 0;
        for (di = 0; di != NUC_NUM_BQS; ++di)
            s->founder_base_marginal[b] += s->complete_jpd[b][di];
    }
    
    for (b = 0; b != 4; ++b)
        for (di = 0; di != NUC_NUM_BQS; ++di)
            s->founder_base_likelihood[b][di] =
                s->complete_jpd[b][di] / s->founder_base_marginal[b];
}


inline int base_to_index(char base)
{
    switch(base)
    {
    case 'a': case 'A': return 0; break;
    case 'c': case 'C': return 1; break;
    case 'g': case 'G': return 2; break;
    case 't': case 'T': return 3; break;
    default: return 4; break;
    }
}


size_t encode_nucleotide(char basecall, size_t quality, size_t strand_index)
{
    assert(quality <= NUC_HIGHEST_QUALITY);
    return 
        (base_to_index(basecall) * NUC_NUM_QS) + (quality * NUC_NUM_S) + strand_index;
}

/* expensive ! */
void decode_nucleotide(size_t code, char * basecall, size_t *quality, size_t *strand)
{
    static char bases_upper[] = "ACGTN";
    *basecall = bases_upper[code / NUC_NUM_QS];
    *quality = (code % NUC_NUM_QS) / NUC_NUM_S; 
    *strand = (code % NUC_NUM_S) == 0 ? NUC_PLUS_STRAND : NUC_MINUS_STRAND;
}


/* initialize fbqs_cpd field of c.  packs the subset of statistics for
   this model that are represented in the c->stats_index */
void nucleotide_stats_pack(const struct nucleotide_stats *stats, struct packed_counts *c)
{
    size_t i, f;
    for (i = 0; i != c->num_data; ++i)
        for (f = 0; f != 4; ++f)
            c->stats[i].cpd[f] = stats->founder_base_likelihood[f][c->stats_index[i]];
}
