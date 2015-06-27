/* */

/* find lowest position p in [beg, end) s.t. elem_is_less(p, par) = 0,
   or end if no such position can be found. */
unsigned virtual_lower_bound(unsigned beg, unsigned end, 
                             int (*less_fcn)(unsigned pos, void *par),
                             void *par)
{
    unsigned mid, len = end - beg;
    while (len > 0)
    {
        mid = beg + (len>>1);
        if (less_fcn(mid, par)) beg = mid + 1;
        else end = mid;
        len = end - beg;
    }
    return beg;
}

/* find the lowest position p in [beg, end) s.t. less_fcn(p, par)
   = 1, or end if no such position can be found */
unsigned virtual_upper_bound(unsigned beg, unsigned end, 
                             int (*less_fcn)(unsigned pos, void *par),
                             void *par)
{
    unsigned mid, len = end - beg;
    while (len > 0)
    {
        mid = beg + (len>>1);
        if (less_fcn(mid, par)) end = mid;
        else beg = mid + 1;
        len = end - beg;
    }
    return beg;
}
