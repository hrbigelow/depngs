#ifndef _PAIR_DIST_STATS_H
#define _PAIR_DIST_STATS_H

#define N_STATS_CATEGORIES 5

#include <unistd.h>

struct pair_dist_stats {
    /* total number of loci pairs for which at least one sample has
       coverage */
    size_t total; 

    /* number of loci in 'total' that fall into each of the categories */
    size_t dist_count[N_STATS_CATEGORIES];

    /* number of dist_count[CHANGED] loci that are further confirmed
       by weighted sampling. */
    size_t confirmed_changed; 

    /* number of loci that are primary or secondary cacheable, based
       on the 8 alpha values.  */
    size_t cacheable;

    /* number of cacheable loci in which the cache was already set. */
    size_t cache_was_set; 
};

/* prints out the internal pair_stats */
void
print_pair_stats(const char *stats_file);


/* increment all pair_dist_stats in bam_sample_pairs with stats */
void
accumulate_pair_stats(struct pair_dist_stats *stats);

#endif /* _PAIR_DIST_STATS_H */
