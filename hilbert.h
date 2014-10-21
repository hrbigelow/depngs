#include <gmp.h>
#include <stdint.h>

void int_to_Hilbert(mpz_t const& index, uint64_t * coords);
void Hilbert_to_int(uint64_t const* coords, mpz_t & index);

void unpack_index(mpz_t const& index, int *chunks);
void pack_index(int const* chunks, mpz_t & index);
void initial_start_end(size_t nchunks, int *start, int *end);
void unpack_coords(const uint64_t *coords, int *unpacked);
void pack_coords(const int *chunks, uint64_t *packed);
int gray_encode(const int &bn);
int gray_decode(const int &gray);
int gray_encode_travel(int start, int end, int mask, int i);
int gray_decode_travel(int start, int end, int mask, int gray);
void child_start_end(int parent_start, int parent_end, int mask, int i,
                     int * child_start, int * child_end);
