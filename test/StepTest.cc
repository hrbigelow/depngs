#include <cstdlib>
#include <cassert>

#include "henry/tools.h"
#include "henry/hilbert.h"
#include "henry/slice_sampling.h"



//Random step test
int main(int argc, char ** argv)
{
    size_t nbits_per_dim = static_cast<size_t>(atoi(argv[1]));
    size_t num_steps = static_cast<size_t>(atoi(argv[2]));
    size_t test_range = static_cast<size_t>(atoi(argv[3]));
    char * initial_point_file = argv[4];

    double initial_point[1000];
    size_t ndim = ParseNumbersFile(initial_point_file, initial_point);

    mpz_t U, N, B, x, xprime, zero;
    mpf_t uniform;

    size_t total_range = nbits_per_dim * ndim;
    
    mpz_init(U);
    mpz_init(N);
    mpz_init(B);
    mpz_init(zero);
    mpz_init(xprime);
    mpf_init(uniform);
    mpz_init(x);

    mpz_set_ui(B, 0);
    mpz_set_ui(zero, 0);
    mpz_set_d(B, pow(2,total_range));

    gmp_randstate_t rand_state;
    gmp_randinit_mt(rand_state);

    uint64_t * xcoord_grid = new uint64_t[ndim];
    double * xcoord = new double[ndim];

    //initialize x
    std::copy(initial_point, initial_point + ndim, xcoord);
    expand_to_grid_coords(xcoord, ndim, xcoord_grid, nbits_per_dim);
    Hilbert_to_int(xcoord_grid, ndim, nbits_per_dim, x);


    for (size_t i = 0; i != num_steps; ++i)
    {
        //goes through the space of integers
        mpz_urandomb(U, rand_state, total_range);
        mpz_urandomb(N, rand_state, test_range);

        mpz_sub(xprime, x, U);

        if (mpz_cmp(xprime, zero) < 0)
        {
            //xprime < 0: wrap it
            mpz_add(xprime, xprime, B);
        }

        mpz_xor(xprime, xprime, N);
        mpz_add(xprime, xprime, U);
        if (mpz_cmp(xprime, B) > 0)
        {
            //xprime > B:  wrap it
            mpz_sub(xprime, xprime, B);
        }

        int_to_Hilbert(xprime, nbits_per_dim, xcoord_grid, ndim);

        contract_from_grid_coords(xcoord_grid, ndim, xcoord, nbits_per_dim);

        if (true)
        {
            printf("%Zu", test_range);
            for (size_t d = 0; d != ndim; ++d)
            {
                printf("\t%20.18f", xcoord[d]);
            }
            printf("\n");
        }

        mpz_set(x, xprime);

    }

    delete xcoord;
    delete xcoord_grid;

    return 0;
}
