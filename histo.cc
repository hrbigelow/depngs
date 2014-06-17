#include "histo.h"

#include <cmath>

// return the bin number that the <val> falls into
// bin 0 is [0, 1/B), center is 1/2B
// bin 1 is [1/B, 2/B), center is 3/2B
// etc
size_t bin(float val, size_t bins_per_unit)
{
    return (size_t) floor(val * bins_per_unit);
}

// return the real coordinate of the center of the bin of a given size
float bin_center(size_t bin, size_t bins_per_unit)
{
    return ((bin * 2l) + 1l) / (2 * (float)bins_per_unit);
}
            

void print_hist(FILE * hist_fh, const char *label, char *contig,
                size_t *hist, size_t bins_per_unit, size_t nbins)
{
    for (size_t bin = 0; bin != nbins; ++bin)
    {
        fprintf(hist_fh, "%s\t%s\t%5.3f\t%Zu\n", label, contig, 
                bin_center(bin, bins_per_unit), hist[bin]);
                        
        hist[bin] = 0;
    }
    fflush(hist_fh);
}
