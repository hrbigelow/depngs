#include "pair_dist_stats.h"

#include <pthread.h>
#include "bam_sample_info.h"
#include "common_tools.h"
#include "binomial_est.h"

static pthread_mutex_t pair_stats_mtx = PTHREAD_MUTEX_INITIALIZER;

/* increment all pair_dist_stats in bam_sample_pairs with stats */
void
accumulate_pair_stats(struct pair_dist_stats *stats)
{
    pthread_mutex_lock(&pair_stats_mtx);
    unsigned s, p;
    for (p = 0; p != bam_sample_pairs.n; ++p)
    {
        bam_sample_pairs.m[p].stats.total += stats[p].total;
        bam_sample_pairs.m[p].stats.confirmed_changed += stats[p].confirmed_changed;
        bam_sample_pairs.m[p].stats.cacheable += stats[p].cacheable;
        bam_sample_pairs.m[p].stats.cache_was_set += stats[p].cache_was_set;

        for (s = 0; s != N_STATS_CATEGORIES; ++s)
            bam_sample_pairs.m[p].stats.dist_count[s] += stats[p].dist_count[s];
    }

    pthread_mutex_unlock(&pair_stats_mtx);
}


void print_pair_stats(const char *stats_file)
{
    FILE *fh = open_if_present(stats_file, "w");
    pthread_mutex_lock(&pair_stats_mtx);
    fprintf(fh, 
            "%s\t%s\t%s\t%s\t%s\t%s", 
            "sample1", "sample2", "total", "cacheable",
            "cache_was_set", "confirmed_changed");

    unsigned s;
    for (s = 0; s != N_STATS_CATEGORIES; ++s)
        fprintf(fh, "\t%s", fuzzy_state_strings[s]);
    fprintf(fh, "\n");

    unsigned p;
    for (p = 0; p != bam_sample_pairs.n; ++p)
    {
        /* Print out statistics */
        int s1 = bam_sample_pairs.m[p].s1,
            s2 = bam_sample_pairs.m[p].s2;
        fprintf(fh, "%s\t%s", 
                bam_samples.m[s1].label,
                s2 == PSEUDO_SAMPLE ? "REF" : bam_samples.m[s2].label);
        fprintf(fh, "\t%zu", bam_sample_pairs.m[p].stats.total);
        fprintf(fh, "\t%zu", bam_sample_pairs.m[p].stats.cacheable);
        fprintf(fh, "\t%zu", bam_sample_pairs.m[p].stats.cache_was_set);
        fprintf(fh, "\t%zu", bam_sample_pairs.m[p].stats.confirmed_changed);

        for (s = 0; s != N_STATS_CATEGORIES; ++s)
            fprintf(fh, "\t%zu", bam_sample_pairs.m[p].stats.dist_count[s]);
        fprintf(fh, "\n");
    }
    fclose(fh);
    pthread_mutex_unlock(&pair_stats_mtx);
}
