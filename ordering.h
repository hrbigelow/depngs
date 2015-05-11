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
int cmp_pair_ordering(const void *pa, const void *pb);

/* ordering of a pair_ordering_range */
int cmp_pair_ordering_range(const void *pa, const void *pb);

#endif /* _ORDERING_H */
