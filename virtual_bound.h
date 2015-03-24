#ifndef _VIRTUAL_BOUND_H
#define _VIRTUAL_BOUND_H

/* find lowest position p in [beg, end) s.t. elem_is_less(p, par) = 0,
   or end if no such position can be found.

   elem_is_less provides both a 'virtual' array of values and a query
   value (through par) to be compared to these virtual values.
   'elem_is_less' returns 1 if the virtual element value at p is less
   than the query value.
  */
unsigned virtual_lower_bound(unsigned beg, unsigned end, 
                             int (*elem_is_less)(unsigned pos, void *par), void *par);

/* find the lowest position p in [beg, end) s.t. query_is_less(p, par)
   = 1 */
unsigned virtual_upper_bound(unsigned beg, unsigned end, 
                             int (*query_is_less)(unsigned pos, void *par), void *par);


#endif /* _VIRTUAL_BOUND_H */
