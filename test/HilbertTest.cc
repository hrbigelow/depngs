#include <cstdlib>
#include <cassert>

#include "henry/hilbert.h"
#include "henry/slice_sampling.h"

//Hilbert space-filling curve test

int main(int argc, char ** argv)
{

    size_t nbits_per_dim = static_cast<size_t>(atof(argv[1]));
    size_t ndim = static_cast<size_t>(atof(argv[2]));

    uint64_t * grid = new uint64_t[ndim];
    double * ucube = new double[ndim];

    uint64_t * grid2 = new uint64_t[ndim];

    size_t i2;
    mpz_t hi, hi2;
    mpz_init(hi);
    mpz_init(hi2);

    int total_bits = nbits_per_dim * ndim;
    for (size_t i = 0; i != static_cast<size_t>(1<<total_bits); ++i)
    {
        mpz_set_ui(hi, i);
        int_to_Hilbert(hi, nbits_per_dim, grid, ndim);
        Hilbert_to_int(grid, ndim, nbits_per_dim, hi2);
        i2 = mpz_get_ui(hi2);

        contract_from_grid_coords(grid, ndim, ucube, nbits_per_dim);

        assert(i == i2);
        printf("%Zu\t%Zu", i, i2);
        for (size_t d = 0; d != ndim; ++d)
        {
            printf("\t%li\t%12.10f", grid[d], ucube[d]);
        }
        printf("\n");
    }

    delete grid;
    delete ucube;
    delete grid2;

    return 0;
}
