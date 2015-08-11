#ifndef _BAM_SAMPLE_INFO_H
#define _BAM_SAMPLE_INFO_H

/* facilities for */

#include <unistd.h>
#include <stdio.h>

/* symbolizes a deeply sequenced sample with genotype matching the
   reference.  used as a 'virtual' sample in the pairwise comparison
   of samples to produce traditional mutation calls. */
#define REFERENCE_SAMPLE UINT_MAX

extern struct bam_sample_info bam_samples;

extern struct bam_sample_pair_info bam_sample_pairs;

/* program-wide initialization of bam sample information.
   sample_pairs_file may be NULL; if so, bam_sample_pairs will be
   initialized with m = NULL, n = 0  */
void bam_sample_info_init(const char *samples_file,
                          const char *sample_pairs_file);

/* call after calling init_sample_pairs and init_sample_attributes */
void bam_sample_info_free();


#define N_STATS_CATEGORIES 5

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


#define PSEUDO_SAMPLE -1


struct bam_sample_pair_info {
    struct { 
        unsigned s1, s2; 
        struct pair_dist_stats stats;
    } *m;
    unsigned n;
};


#define MAX_LABEL_LEN 100

struct bam_sample_info {
    struct {
        char *bam_file;
        char label[MAX_LABEL_LEN + 1];
    } *m;
    unsigned n;
};


/* prints out the internal pair_stats */
void
print_pair_stats(const char *stats_file);


/* increment all pair_dist_stats in bam_sample_pairs with stats */
void
accumulate_pair_stats(struct pair_dist_stats *stats);


#endif /* _BAM_SAMPLE_INFO_H */
