/* 
   Routines to identify the G nearest neighbors of a point p in N
dimensions.

 */
#include <stdlib.h>
#include <float.h>
#include <math.h>
#include <assert.h>
#include <stdio.h>
#include "spatial_search.h"

static unsigned cmp_dim;

void set_cmp_dim(unsigned c)
{
    cmp_dim = c;
}

/* comparison by marked_point::p[cmp_dim].  cmp_dim is external */
int marked_point_comp(const void *pa, const void *pb)
{
    const struct marked_point
        *a = *(struct marked_point **)pa,
        *b = *(struct marked_point **)pb;

    return
        a->p[cmp_dim] < b->p[cmp_dim] ? -1 : a->p[cmp_dim] > b->p[cmp_dim] ? 1 : 0;
}


double euclidean_dist(double *p1, double *p2)
{
    int d;
    double sq_dist = 0;

    for (d = 0; d != NDIM; ++d)
        sq_dist += pow(p1[d] - p2[d], 2.0);

    return sqrt(sq_dist);
}


/* returns the left-most pointer p in [base, base + nmemb * size)
   where compar(p, val) >= 0 */
void *lower_bound(void *base, size_t nmemb, size_t size, 
                  void *val, /* address of the value to be fed into the comparison function */
                  int(*compar)(const void *, const void *))
{
    size_t half;
    char *first = (char *)base, *middle;
    while (nmemb > 0)
    {
        half = nmemb >> 1;
        middle = first;
        middle += half * size;
        if (compar(middle, val) < 0)
        {
            first = middle;
            first += size;
            nmemb = nmemb - half - 1;
        }
        else
            nmemb = half;
    }
    return first;
}


/* returns the right-most pointer p in [base, base + nmemb * size)
   where compar(val, p) < 0 or p == base + nmemb * size if this is not
   true of any element. */
void *upper_bound(void *base, size_t nmemb, size_t size, 
                  void *val, /* address of the value to be fed into the comparison function */
                  int(*compar)(const void *, const void *))
{
    size_t half;
    char *first = (char *)base, *middle;
    while (nmemb > 0)
    {
        half = nmemb >> 1;
        middle = first;
        middle += half * size;
        if (compar(val, middle) < 0)
            nmemb = half;
        else
        {
            first = middle;
            first += size;
            nmemb = nmemb - half - 1;
        }
    }
    return first;
}




/* among the left and right points, find the closest point in the
   single-component distance of 'dim' dimension to 'center'.  set
   '*closest_point' equal to this found point, and return its linear
   distance to 'center'
  */
double min_dist(struct marked_point **left[],
                struct marked_point **right[],
                struct marked_point *center,
                struct marked_point **closest_point,
                int *min_dim)
{
    int d;
    double min_dist = DBL_MAX, cur_dist;
    for (d = 0; d != NDIM; ++d)
    {   
        if ((cur_dist = center->p[d] - (*left[d])->p[d]) < min_dist)
        {
            min_dist = cur_dist;
            *min_dim = d;
            *closest_point = *left[d];
        }
        if ((cur_dist = (*right[d])->p[d] - center->p[d]) < min_dist)
        {
            min_dist = cur_dist;
            *min_dim = d;
            *closest_point = *right[d];
        }
    }
    return min_dist;
}


int all_flags_set(unsigned char *flags)
{
    int i, all_set = 1;
    for (i = 0; i != NDIM; ++i)
        all_set &= flags[i];
    return all_set;
}


#define FOUR_THIRDS_PI (M_PI * 4.0 / 3.0)

/* estimate the volume associated with a particular center
   point. 'head' is the head of a linked list of the G nearest
   neighbors that terminates at the center point in question */
double estimate_volume(struct marked_point *head, unsigned G)
{
    /* struct marked_point *p = head; */
    /* unsigned g = G; */
    /* while (p->radius) */
    /* { */
    /*     double est_radius = (p->radius + p->prev->radius) / 2.0; */
    /*     fprintf(stdout, "%i\t%20.18f\n", g, */
    /*             FOUR_THIRDS_PI * pow(est_radius, NDIM) / (double)g); */
    /*     p = p->prev; */
    /*     --g; */
    /* } */
    /* fprintf(stdout, "\n\n"); */
    double est_radius = (head->radius + head->prev->radius) / 2.0;
    return FOUR_THIRDS_PI * pow(est_radius, NDIM) / (double)(G - 1);
}


/* returns the head of a linked list of the G nearest neighbors of a
   given point. it is the caller's responsibility to clear the in_grid
   flags of all points */
struct marked_point *spatial_search(struct marked_point **points[NDIM],
                                    unsigned npoints,
                                    struct marked_point *center,
                                    int G)
{
    double grid_width = 0;
    int num_nbor = 1; /* we count the head_node as a neighbor initially */
    struct marked_point head_node;
    head_node.radius = DBL_MAX;
    head_node.prev = center;

    center->radius = 0;
    
    struct marked_point *cur, *head = &head_node;
    struct marked_point **left[NDIM], **right[NDIM], **pcenter;

    /* initialize left and right */
    int d;
    for (d = 0; d != NDIM; ++d)
    {
        set_cmp_dim(d);
        pcenter = (struct marked_point **)lower_bound(points[d], npoints, 
                                                      sizeof(struct marked_point *),
                                                      &center, marked_point_comp);
        while (*pcenter != center)
            ++pcenter;

        left[d] = pcenter == points[d] ? pcenter : pcenter - 1;
        right[d] = pcenter == points[d] + npoints - 1 ? pcenter : pcenter + 1;
    }

    int min_dim;
    
    while (grid_width < head->radius)
    {
        grid_width = min_dist(left, right, center, &cur, &min_dim);

        /* update left or right */
        if (cur == *left[min_dim])
        {
            if (left[min_dim] == points[min_dim])
                return NULL;
            --left[min_dim];
        }
        else
        {
            if (right[min_dim] == points[min_dim] + npoints - 1)
                return NULL;
            ++right[min_dim];
        }

        if (cur->center != center)
        {
            cur->center = center;
            for (d = 0; d != NDIM; ++d)
                cur->in_grid[d] = 0;
            cur->radius = 0;
        }
        cur->in_grid[min_dim] = 1;

        if (all_flags_set(cur->in_grid))
            cur->radius = euclidean_dist(cur->center->p, cur->p);

        if (cur->radius && cur->radius < head->radius)
        {
            struct marked_point *p = head;
            while (cur->radius < p->prev->radius)
                p = p->prev;

            /* splice */
            cur->prev = p->prev;
            p->prev = cur;
            ++num_nbor;
            assert(cur->prev != cur);

            /* shrink the list */
            if (num_nbor > G)
            {
                head = head->prev;
                --num_nbor;
            }
        }
    }
    return head;
}
