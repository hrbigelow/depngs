#ifndef _JOINT_BOUND_SEARCH_H
#define _JOINT_BOUND_SEARCH_H

#include "ordering.h"

typedef size_t (*range_to_size_func_t)(void *index, 
                                       struct pair_ordering beg,
                                       struct pair_ordering end);

typedef struct pair_ordering
(*size_to_pos_func_t)(void *index, 
                      struct pair_ordering beg,
                      size_t size);

/* the main configuration object. opaque type using pimpl pattern */
struct joint_bound_search;

/* initialize the main configuration object */
struct joint_bound_search *
joint_bound_search_init(pos_to_off_func_t _size_func,
                        pos_to_size_func_t _range_func,
                        void *indexes,
                        size_t num_files,
                        size_t mem);

/* free all resources */
void joint_bound_search_free(struct joint_bound_search *jb);

/* returns the position after start_pos closest to using a total of
   jb->mem.  sets lb_off and ub_off with the offset ranges
   corresponding with start_pos and the returned upper bound
   position. */
struct pair_ordering
joint_bound_search(struct joint_bound_search *jb,
                   struct pair_ordering pos_start);

#endif /* _JOINT_BOUND_SEARCH_H */
