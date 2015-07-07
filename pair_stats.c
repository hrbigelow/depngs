#include <pthread.h>


static pthread_mutex_t pair_stats_mtx = PTHREAD_MUTEX_INITIALIZER;

#define N_STATS_CATEGORIES \
    sizeof(sample_pairs.p[0].stats.dist_count) \
    / sizeof(sample_pairs.p[0].stats.dist_count[0])


/* increment all pair_dist_stats in sample_pairs with stats */
void accumulate_pair_stats(struct pair_dist_stats *stats)
{
    pthread_mutex_lock(&pair_stats_mtx);
    unsigned s, p;
    for (p = 0; p != sample_pairs.n; ++p)
    {
        sample_pairs.p[p].stats.total += stats[p].total;
        sample_pairs.p[p].stats.confirmed_changed += stats[p].confirmed_changed;
        sample_pairs.p[p].stats.cacheable += stats[p].cacheable;
        sample_pairs.p[p].stats.cache_was_set += stats[p].cache_was_set;

        for (s = 0; s != N_STATS_CATEGORIES; ++s)
            sample_pairs.p[p].stats.dist_count[s] += stats[p].dist_count[s];
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
    for (p = 0; p != sample_pairs.n; ++p)
    {
        /* Print out statistics */
        int s1 = sample_pairs.p[p].s1,
            s2 = sample_pairs.p[p].s2;
        fprintf(fh, "%s\t%s", 
                samples.atts[s1].label,
                s2 == PSEUDO_SAMPLE ? "REF" : samples.atts[s2].label);
        fprintf(fh, "\t%zu", sample_pairs.p[p].stats.total);
        fprintf(fh, "\t%zu", sample_pairs.p[p].stats.cacheable);
        fprintf(fh, "\t%zu", sample_pairs.p[p].stats.cache_was_set);
        fprintf(fh, "\t%zu", sample_pairs.p[p].stats.confirmed_changed);

        for (s = 0; s != N_STATS_CATEGORIES; ++s)
            fprintf(fh, "\t%zu", sample_pairs.p[p].stats.dist_count[s]);
        fprintf(fh, "\n");
    }
    fclose(fh);
    pthread_mutex_unlock(&pair_stats_mtx);
}
