#include <immintrin.h>

#if 0
/* write a vectorized version of frexp, except that only handles normals */
double vec_frexp (double x, int *e)
{
    int32_t hx, ix, lx;

    /* PSHUFD (Shuffle packed double-words) (perform twice) */
    EXTRACT_WORDS (hx, lx, x);

    /* VANDPS (bitwise logical and.  not clear why this would care about integers... */
    ix = 0x7fffffff & hx;

    /* PSRLD (shift each packed double-word right by n bits */
    /* VPSUBD (subtract packed double-word integers in-place) */
    *e = (ix >> 20) - 1022;     

    /* VANDPS (zero-out the exponent) */
    /* VPOR (set exponent to 1022) */
    hx = (hx & 0x800fffff) | 0x3fe00000;

    /* PSHUFD (Shuffle packed double-words) */
    SET_HIGH_WORD (x, hx);

    return x;
}
#endif

#ifdef __AVX2__
void vec_frexp(__m256d x1, __m256d x2, __m256d *n1, __m256d *_n2, __m256i *e)
{
    __m256i twiddle = _mm256_set_epi32(1, 0, 3, 2, 5, 4, 7, 6);
    __m256i lohi = _mm256_permutevar8x32_epi32((__m256i)x2, twiddle); /* LHLHLHLH */
    __m256i lo = _mm256_blend_epi32(x1, lohi, 0x99);
    __m256i hi = _mm256_blend_epi32(lohi, x1, 0x99);
    __m256i exp = _mm256_and_si256(hi, _mm256_set1_epi32(0x7fffffff));
    exp = _mm256_srli_si256(exp, 20);
    exp = _mm256_sub_epi32(exp, _mm256_set1_epi32(1022));
    __m256i hi2 = _mm256_and_si256(hi, _mm256_set1_epi32(0x800fffff));
    hi2 = _mm256_or_si256(hi2, _mm256_set1_epi32(0x3fe00000));
    *n1 = _mm256_blend_epi32(hi2, lo, 0x99);
    __m256i trg2 = _mm256_blend_epi32(lo, hi2, 0x99);
    *n2 = _mm256_permutevar8x32_epi32(trg2, twiddle);
}
#endif

#ifdef __AVX__
void vec_frexp(__m256d x1, __m256d x2, __m256d *n1, __m256d *_n2, __m256i *e)
{
    __m256i twiddle = _mm256_set_epi32(1, 0, 3, 2, 5, 4, 7, 6);
    __m256i lohi = (__m256i)_mm256_permutevar_ps((__m256)x2, twiddle); /* LHLHLHLH */
    __m256i lo = (__m256i)_mm256_blend_ps(x1, lohi, 0x55); /* combine HLHLHLHL and LHLHLHLH, using 01010101*/
    __m256i hi = (__m256i)_mm256_blend_ps(lohi, x1, 0x55);
    __m256i exp = (__m256i)_mm256_and_pd((__m256d)hi, (__m256d)_mm256_set1_epi32(0x7fffffff));
    exp = _mm256_srli_si256(exp, 20);
    exp = _mm256_sub_epi32(exp, _mm256_set1_epi32(1022));
    __m256i hi2 = _mm256_and_si256(hi, _mm256_set1_epi32(0x800fffff));
    hi2 = _mm256_or_si256(hi2, _mm256_set1_epi32(0x3fe00000));
    *n1 = _mm256_blend_epi32(hi2, lo, 0x99);
    __m256i trg2 = _mm256_blend_epi32(lo, hi2, 0x99);
    *n2 = _mm256_permutevar8x32_epi32(trg2, twiddle);
}
#endif


/* 
src1 HLHLHLHL
src2 LHLHLHLH
imm  10101010
lo   LLLLLLLL
hi   HHHHHHHH
exp  eeeeeeee
hi2  hhhhhhhh  -- the normalized high double-word
trg1 hLhLhLhL  -- blend
trg2 LhLhLhLh  -- blend
hilo hLhLhLhL  -- permute


PSRLDQ shifts by number of bytes.
SHRD multiprecision shifts of 64 bits or more.
VPSRLVD  (AVX2)

PALIGNR
PBLENDW
PEXTRW
PMOVSX
POR
PAND
PSHUFHW
PSUB
SAL/SAR (shift)
SHLD
SHUFPD
VEXTRACTI128
VFMADD132SD
VGATHERDPD
VPBLENDD
VPERMPD
VZEROUPPER

*/
