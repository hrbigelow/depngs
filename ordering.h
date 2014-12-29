#ifndef _ORDERING_H
#define _ORDERING_H

#include <stdint.h>
#include <unistd.h>

/* a generic space in which to store line ordinal information */
struct pair_ordering {
    size_t hi, lo;
};

struct pair_ordering_range {
    struct pair_ordering beg, end;
};

/* defines the ordering logic for the pair_ordering structure. */
int less_pair_ordering(const void *pa, const void *pb);

/* ordering of a pair_ordering_range */
int less_pair_ordering_range(const void *pa, const void *pb);


static const struct pair_ordering min_pair_ord = { 0, 0 };
static const struct pair_ordering max_pair_ord = { SIZE_MAX, SIZE_MAX };

#endif /* _ORDERING_H */
