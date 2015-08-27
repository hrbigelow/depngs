#include "ordering.h"

#define CMP(a, b) ((a) < (b) ? -1 : (a) > (b) ? 1 : 0)

/* defines the ordering logic for the pair_ordering structure. */
int cmp_pair_ordering(const void *pa, const void *pb)
{
    const struct pair_ordering *a = pa, *b = pb;

    int acmp;
    return
        (acmp = CMP(a->hi, b->hi)) != 0
        ? acmp
        : CMP(a->lo, b->lo);
}

/* ordering of a pair_ordering_range.   */
int cmp_pair_ordering_range(const void *pa, const void *pb)
{
    const struct pair_ordering_range *a = pa, *b = pb;

    int bcmp;
    return 
        (bcmp = cmp_pair_ordering(&a->beg, &b->beg)) != 0
        ? bcmp
        : cmp_pair_ordering(&a->end, &b->end);
}


/* compare two ranges in the following way, based on the comparison of
   the beg, end of each range:
   0  0:  0
   0  1:  0
   0 -1:  0
   1  0:  0
  -1  0:  0
   1 -1:  0
  -1  1:  0
   1  1:  1
  -1 -1: -1
*/
int cmp_pair_ordering_contained(const void *pa, const void *pb)
{
    const struct pair_ordering_range *a = pa, *b = pb;
    int bc = cmp_pair_ordering(&a->beg, &b->beg);
    int ec = cmp_pair_ordering(&a->end, &b->end);
    if (bc == ec && bc != 0) return bc;
    else return 0;
}
