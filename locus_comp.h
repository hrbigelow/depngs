#ifndef _LOCUS_COMP_H
#define _LOCUS_COMP_H

#include <cstring>
#include <map>

struct ltstr
{
    bool operator()(const char* s1, const char* s2) const
    {
        return strcmp(s1, s2) < 0;
    }
};

struct less_locus_position
{
    std::map<char const*, size_t, ltstr> * contig_order;
    bool operator()(char const* locus_line1, char const* locus_line2) const;
};

struct equal_locus_position
{
    std::map<char const*, size_t, ltstr> * contig_order;
    bool operator()(char const* locus_line1, char const* locus_line2) const;
};


#endif // _LOCUS_COMP_H
