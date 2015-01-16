#include "nucleotide_stats.h"

#include <stdio.h>
#include <string.h>
#include <assert.h>

void nucleotide_stats_initialize(const char *rdb_file, struct nucleotide_stats *s)
{

    FILE *rdb_fh = fopen(rdb_file, "r");
    if (rdb_fh == NULL)
    {
        fprintf(stderr, "Couldn't open data jpd file %s\n", rdb_file);
        exit(5);
    }

    double ct[4];
    char base, strand;
    int qual;
    size_t bqs, f;

    memset(s->marg, 0, sizeof(s->marg));

    while (! feof(rdb_fh))
    {
        (void)fscanf(rdb_fh, "%c_%i_%c\t%lf\t%lf\t%lf\t%lf\n", 
                     &base, &qual, &strand, ct, ct + 1, ct + 2, ct + 3);
        
        if (ct[0] < 0 || ct[1] < 0 || ct[2] < 0 || ct[3] < 0)
        {
            fprintf(stderr, "%s: found negative count for %c_%i_%c.\n", 
                    __func__, base, qual, strand);
            exit(11);
        }
        if (ct[0] + ct[1] + ct[2] + ct[3] == 0) continue;
        bqs = 
            encode_nucleotide(base, qual, 
                              (strand == '+' ? NUC_PLUS_STRAND : NUC_MINUS_STRAND));
        
        for (f = 0; f != 4; ++f)
        {
            s->jpd[f][bqs] = ct[f];
            s->cpd[f][bqs] = ct[f];
            s->marg[f] += ct[f];
        }
    }
    fclose(rdb_fh);

    /* normalize jpd_buffer */
    double sum_inv = 1.0 / (s->marg[0] + s->marg[1] + s->marg[2] + s->marg[3]);

    double marg_inv[] = 
        { 1.0 / s->marg[0], 1.0 / s->marg[1], 1.0 / s->marg[2], 1.0 / s->marg[3] };

    /* normalize marg */
    for (f = 0; f != 4; ++f) s->marg[f] *= sum_inv;

    for (f = 0; f != 4; ++f)
        for (bqs = 0; bqs != NUC_NUM_BQS; ++bqs)
        {
            s->jpd[f][bqs] *= sum_inv;
            s->cpd[f][bqs] *= marg_inv[f];
        }
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
void nucleotide_stats_pack(const struct nucleotide_stats *stats, struct packed_counts *pc)
{
    size_t i, f;
    for (i = 0; i != pc->num_data; ++i)
        for (f = 0; f != 4; ++f)
            pc->stats[i].cpd[f] = stats->cpd[f][pc->stats_index[i]];
}
