#ifndef _FOUNDER_DIST_H
#define _FOUNDER_DIST_H

#define MAX_QUAL_BINS 9

struct founder_pair_counts {
    uint16_t major, minor;
};

struct binned_qual_counts {
    struct founder_pair_counts bins[MAX_QUAL_BINS];
};

#endif /* _FOUNDER_DIST_H */
