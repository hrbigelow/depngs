#ifndef _SPATIAL_SEARCH_H
#define _SPATIAL_SEARCH_H

#define NDIM 4

struct marked_point {
    double p[NDIM]; /* the underlying spatial coordinates of this point */
    double dist_val; /* value of the underlying distribution at this point */
    unsigned char in_grid[NDIM]; /* true if this point is in the grid
                                    at that dimension */
    struct marked_point *center; /* the center point for this search */
    struct marked_point *prev; /* points to the previous marked_point
                                  in the linked list of up to G
                                  nearest neighbors */
    double radius; /* the distance from this point to the originating
                      point */
};


void set_cmp_dim(unsigned c);

/* comparison by marked_point::p[cmp_dim].  cmp_dim is external */
int marked_point_comp(const void *pa, const void *pb);


/* estimate the volume associated with a particular center
   point. 'head' is the head of a linked list of the G nearest
   neighbors that terminates at the center point in question */
double estimate_volume(struct marked_point *head, unsigned G);

/* returns the head of a linked list of the G nearest neighbors of a
   given point. it is the caller's responsibility to clear the in_grid
   flags of all points.  return NULL if the search hit a boundary. */
struct marked_point *
spatial_search(struct marked_point **points[NDIM],
               unsigned npoints,
               struct marked_point *center,
               int G);


#endif /* _SPATIAL_SEARCH_H */
