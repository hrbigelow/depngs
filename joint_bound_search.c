#include "joint_bound_search.h"

#include <stdlib.h>

/* Traverse N files together incrementally.  At each increment, find a
   single upper-bound position, and the physical offset in each file
   corresponding to that upper-bound position, such that the total
   size of physical offset differences is close to a desired amount
   'mem'.
*/

struct joint_bound_search
{
    range_to_size_func_t find_size;
    size_to_pos_func_t find_range;
    void *indexes; /* set of 'num_files' index structures used for searching */
    size_t num_files;
    size_t mem;
};


/* initialize the main configuration object */
struct joint_bound_search *
joint_bound_search_init(range_to_size_func_t _find_size,
                        size_to_pos_func_t _find_range,
                        void *indexes,
                        size_t num_files,
                        size_t mem)
{
    struct joint_bound_search *jb = malloc(sizeof(struct joint_bound_search));
    jb->find_size = _find_size;
    jb->find_range = _find_range;
    jb->indexes = indexes;
    jb->num_files = num_files;
    jb->mem = mem;
    return jb;
}

/* free all resources */
void joint_bound_search_free(struct joint_bound_search *jb)
{
    free(jb);
}


/* this does the expected intercalation when beg.hi == end.hi.  When
   they are different, we have the problem that we don't know the
   maximal value of beg.lo, and so cannot know what the actual point
   of intercalation is.  In this case, we use a two step approach.
   First, guess that the point is somewhere on end.hi, using end.lo =
   0 as a signal.  If that is not true, jump forward on beg.hi by a
   modest fraction until we overshoot. */
struct pair_ordering
intercalate(struct pair_ordering beg, struct pair_ordering end, float frac)
{
    /* meant to grow reasonably fast. */
    static const float growth_factor = 1.1;
    struct pair_ordering mid;
    if (beg.hi == end.hi)
    {
        mid.hi = beg.hi;
        mid.lo = beg.lo + (end.lo - beg.lo) * frac;
    }
    else
    {
        if (end.lo == 0)
        {
            mid.hi = end.hi;
            mid.lo = 0;
        }
        else
        {
            mid.hi = beg.hi;
            mid.lo = beg.hi * growth_factor;
        }
    }

    return mid;
}
            

/* Do the search.  Proceeds as follows:

   1. Speculatively find the N end positions P[i] (one for each file)
      that fall within mem/N size.

   2. It then takes the min, mid, and max among these end positions
      and calculates the sum of sizes of ranges [pos_start, min) etc
      for each file.
      
   3. Start the binary search.  Depending on whether the middle
      sum-of-sums (sum_mid) is greater or less than mem, shrink the
      search interval to [min, mid) or [mid, max).

   4. Find a midpoint position mid using intercalation.

   5. Calculate sum_mid over all files of the range [pos_start, mid),
      and repeat from step 3 until convergence.
 */
struct pair_ordering
joint_bound_search(struct joint_bound_search *jb,
                   struct pair_ordering pos_start)
{

    struct pair_ordering *pos_end = 
        malloc(jb->num_files * sizeof(struct pair_ordering));

    
    size_t 
    void **index = (void **)jb->indexes;
    size_t i;
    for (i = 0; i != jb->num_files; ++i)
        pos_end[i] = jb->find_range(index[i], pos_start, jb->mem / jb->num_files);
    
    qsort(pos_end, jb->num_files, sizeof(struct pair_ordering), less_pair_ordering);

    size_t sum_max = 0, sum_min = 0, sum_mid = 0;

    struct pair_ordering 
        pos_min = pos_end[0],
        pos_mid = pos_end[jb->num_files / 2],
        pos_max = pos_end[jb->num_files - 1];

    for (i = 0; i != jb->num_files; ++i)
        sum_min += jb->find_size(index[i], pos_start, pos_min);

    if (less_pair_ordering(pos_min, pos_mid) != 0)
        for (i = 0; i != jb->num_files; ++i)
            sum_mid += jb->find_size(index[i], pos_start, pos_mid);
    else sum_mid = sum_min;

    if (less_pair_ordering(pos_mid, pos_max) != 0)
        for (i = 0; i != jb->num_files; ++i)
            sum_max += jb->find_size(index[i], pos_start, pos_max);
    else sum_max = sum_mid;

    /* Iterate until convergence. */
    float frac;
    while (sum_min != sum_max)
    {
        /* shrink interval bounds */
        if (sum_mid < jb->mem) 
            sum_min = sum_mid, pos_min = pos_mid;
        else 
            sum_max = sum_mid, pos_max = pos_mid;

        frac = (float)(jb->mem - sum_min) / (float)(sum_max - sum_min);
        pos_mid = intercalate(pos_min, pos_max, frac);

        /* compute new sum */
        sum_mid = 0;
        for (i = 0; i != jb->num_files; ++i)
            sum_mid += jb->find_size(index[i], pos_start, pos_mid);
    }
    free(pos_end);
    return pos_mid;
}
