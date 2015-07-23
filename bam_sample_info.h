#ifndef _BAM_SAMPLE_INFO_H
#define _BAM_SAMPLE_INFO_H

#include <unistd.h>
#include <stdio.h>

#include "pair_dist_stats.h"

extern struct bam_sample_info bam_samples;

extern struct bam_sample_pair_info bam_sample_pairs;

/* call before calling init_sample_pairs and init_sample_attributes */
void bam_sample_info_init(const char *samples_file,
                          const char *sample_pairs_file);

/* call after calling init_sample_pairs and init_sample_attributes */
void bam_sample_info_free();


#define PSEUDO_SAMPLE -1


struct bam_sample_pair_info {
    struct { 
        int s1, s2; 
        struct pair_dist_stats stats;
    } *m;
    unsigned n;
};


#define MAX_LABEL_LEN 100

struct bam_sample_info {
    struct {
        char *bam_file;
        char label[MAX_LABEL_LEN + 1];
        FILE *fh;
    } *m;
    unsigned n;
};



#endif /* _BAM_SAMPLE_INFO_H */
