#include "locus_comp.h"

#include <cstdio>
#include <cassert>

bool less_locus_position::operator()(char const* locus_line1, char const* locus_line2) const
{
    char contig1[100];
    char contig2[100];
    size_t position1;
    size_t position2;
    sscanf(locus_line1, "%s\t%zu", contig1, &position1);
    sscanf(locus_line2, "%s\t%zu", contig2, &position2);

    if (strcmp(contig1, contig2) == 0) 
    { 
        return position1 < position2; 
    }
    else
    {
        std::map<char const*, size_t>::iterator it1, it2;
        it1 = this->contig_order->find(contig1);
        it2 = this->contig_order->find(contig2);
        assert(it1 != this->contig_order->end());
        assert(it2 != this->contig_order->end());
        return (*it1).second < (*it2).second;
    }
}



bool equal_locus_position::operator()(char const* locus_line1, char const* locus_line2) const
{
    char contig1[100];
    char contig2[100];
    size_t position1;
    size_t position2;
    sscanf(locus_line1, "%s\t%zu", contig1, &position1);
    sscanf(locus_line2, "%s\t%zu", contig2, &position2);
    return strcmp(contig1, contig2) == 0 && position1 == position2;
}
