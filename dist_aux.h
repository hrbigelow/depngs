#ifndef _DIST_AUX_H
#define _DIST_AUX_H

#include <unistd.h>
#include <stdio.h>

#include "nucleotide_stats.h"

struct pair_dist_stats {
    size_t total; /* total number of loci pairs for which at least one
                     sample has coverage */
    size_t dist_count[5]; /* number of loci in 'total' that fall into
                             each of the categories */
    size_t confirmed_changed; /* number of dist_count[CHANGED] loci
                                 that are further confirmed by
                                 weighted sampling. */
    size_t cacheable; /* number of loci that are primary or secondary
                         cacheable, based on the 8 alpha values.  */
    size_t cache_was_set; /* number of cacheable loci in which the
                             cache was already set. */
};

/* call before calling init_sample_pairs and init_sample_attributes */
void init_dist_aux();

/* call after calling init_sample_pairs and init_sample_attributes */
void free_dist_aux();

struct sample_pairs {
    struct { 
        int s1, s2; 
        struct pair_dist_stats stats;
    } *p;
    unsigned n;
};

struct sample_pairs
init_sample_pairs(const char *pair_file);

#define MAX_LABEL_LEN 100


struct sample_attrs {
    struct {
        char *file;
        char label[MAX_LABEL_LEN + 1];
        struct nucleotide_stats nuc_stats;
        FILE *fh;
    } *atts;
    unsigned n;
};


struct sample_attrs
init_sample_attributes(const char *samples_file);



#endif /* _DIST_AUX_H */
