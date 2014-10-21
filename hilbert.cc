/*
  implementation of http://www.tiac.net/~sw/2008/10/Hilbert/hilbert.py

  From the website:
  
  To encode from index to coordinates:
  
  * "Unpack" the index into a list of h D-bit ints (called "index
  * chunks").
  
  * h will be the number of levels of recursion (in our case the
  * number of times to loop); from this calculate the orientation of
  * the overall cube.
  
  Then, for each index chunk, starting at the most- significant,
  
  * Use a modified Gray code to convert the index chunk into a D-bit
  * "coordinate chunk" with one bit destined for each of the D
  * coordinates;
  
  * Calculate the start and end corners for the child cube to which
  * the indexed point belongs;
  
  * Loop to do the child cube.
  
  Finally,
  
  * "Pack" the output coordinates by redistributing the D bits from
  * each of h coordinate chunks into D ints, each with h bits.

*/

#include "hilbert.h"
#include "slice_sampling.h"

#include <algorithm>
#include <gmp.h>

#define NCHUNKS NBITS_PER_DIM

void unpack_index(mpz_t const& index, int *chunks);
void pack_index(const int *chunks, mpz_t & index);
void initial_start_end(int *start, int *end);
void unpack_coords(const uint64_t *coords, int *unpacked);
void pack_coords(int const* chunks, uint64_t *packed);
int gray_encode(int const& bn);
int gray_decode(int const& gray);
int gray_encode_travel(int start, int end, int mask, int i);
int gray_decode_travel(int start, int end, int mask, int gray);
void child_start_end(int parent_start, int parent_end, int mask, int i,
                     int *child_start, int *child_end);


//map index to the NDIM-dimensional uint64_t hypercube hilbert coordinate
void int_to_Hilbert(mpz_t const& index, uint64_t * coords)
{
    int index_chunks[NCHUNKS];
    int coord_chunks[NCHUNKS];

    unpack_index(index, index_chunks);
    int mask = (1<<NDIM) - 1;
    int start, end;
    initial_start_end(&start, &end);
    for (size_t j = 0; j != NCHUNKS; ++j)
    {
        int i = index_chunks[j];
        coord_chunks[j] = gray_encode_travel(start, end, mask, i);
        child_start_end(start, end, mask, i, &start, &end);
    }
    pack_coords(coord_chunks, coords);

}


//map NDIM-dimensional uint64_t hypercube hilbert coordinate to (64*NDIM) - bit index on integer line
void Hilbert_to_int(uint64_t const* coords, mpz_t & index)
{
    int index_chunks[NCHUNKS];
    int coord_chunks[NCHUNKS];
    
    unpack_coords(coords, coord_chunks);

    int mask = (1<<NDIM) - 1;
    int start, end;
    initial_start_end(&start, &end);
    for (size_t j = 0; j != NCHUNKS; ++j)
    {
        int i = gray_decode_travel(start, end, mask, coord_chunks[j]);
        index_chunks[j] = i;
        child_start_end(start, end, mask, i, &start, &end);
    }
    pack_index(index_chunks, index);
}



void initial_start_end(int *start, int *end)
{
    *start = 0;
    *end = 1<<(abs(- (NCHUNKS) - 1) % NDIM);
}


void unpack_index(mpz_t const& index, int * chunks)
{
    mpz_t i, zchunk, p;

    mpz_init_set(i, index);
    mpz_init(zchunk);
    mpz_init_set_ui(p, 1<<NDIM);

    for (size_t j = NCHUNKS; j != 0; --j)
    {
        mpz_mod(zchunk, i, p);
        chunks[j-1] = mpz_get_ui(zchunk);
        mpz_div(i, i, p);
    }        
}


void pack_index(const int *chunks, mpz_t & index)
{
    mpz_t p;
    mpz_init(index);
    mpz_init_set_ui(p, 1<<NDIM);

    for (size_t j = 0; j != NCHUNKS; ++j)
    {
        mpz_mul(index, index, p);
        mpz_add_ui(index, index, chunks[j]);
    }

}


//lines up each coord as a row from top to bottom
//loads 'unpacked' with the columns read top-to-bottom, in order from left to right
void unpack_coords(uint64_t const* coords, int * unpacked)
{
    std::fill(unpacked, unpacked + NCHUNKS, 0);
    
    size_t rrev;

    for (size_t r = 0; r != NDIM; ++r)
    {
        rrev = NDIM - r - 1;
        for (size_t c = 0; c != NCHUNKS; ++c)
        {
            unpacked[NCHUNKS - c - 1] |= (static_cast<uint64_t>(((coords[r] >> c) & 1)) << rrev);
        }
    }
}


void pack_coords(int const* chunks, uint64_t * packed)
{
    std::fill(packed, packed + NDIM, 0);
    
    size_t rrev;
    
    for (size_t r = 0; r != NCHUNKS; ++r)
    {
        rrev = NCHUNKS - r - 1;
        for (size_t c = 0; c != NDIM; ++c)
        {
            packed[NDIM - c - 1] |= (((uint64_t)((chunks[r] >> c) & 1)) << rrev);
        }
    }
}


int gray_encode(int const& bn)
{
    return bn ^ (bn / 2);
}


int gray_decode(int const& n)
{
    char div;
    int bn = n, sh = 1;
    while (1)
    {
        div = bn >> sh;
        bn ^= div;
        if (div <= 1)
        {
            return bn;
        }
        sh <<= 1;
    }
}


int gray_encode_travel(int start, int end, int mask, int i)
{
    int travel_bit = start ^ end;
    int modulus = mask + 1;
    int gray = gray_encode(i) * (travel_bit * 2);
    return ((gray | (gray / modulus)) & mask) ^ start;
}


int gray_decode_travel(int start, int end, int mask, int gray)
{
    int travel_bit = start ^ end;
    int modulus = mask + 1;
    int rgray = (gray ^ start) * (modulus / (travel_bit * 2));
    return gray_decode((rgray | (rgray / modulus)) & mask);
}



void child_start_end(int parent_start, int parent_end, int mask, int i,
                     int * child_start, int * child_end)
{
    int start_i = std::max(0, (i-1) & ~1);
    int end_i = std::min(mask, (i+1) | 1);
    *child_start = gray_encode_travel(parent_start, parent_end, mask, start_i);
    *child_end = gray_encode_travel(parent_start, parent_end, mask, end_i);
}
