/* A slightly more efficient implementation of gsl_pow_int */
#include <immintrin.h>

double pow_uint(double x, unsigned int n)
{
    double v = x;
    n >>= 1;
    while (n)
    {
        x *= x;
        if (n & 1) v *= x;
        n >>= 1;
    }
    return v;
}



#ifdef __AVX__

__m256d pow_uint_vec(__m256d x, __m256i n)
{
    __m128i hn, ln; /* high quadword of n */
    ln = _mm_srli_si128(_mm256_extractf128_si256(n, 0), 1);
    hn = _mm_srli_si128(_mm256_extractf128_si256(n, 1), 1);
    n = _mm256_insertf128_si256(n, ln, 0);
    n = _mm256_insertf128_si256(n, hn, 1);
    __m256i zero = _mm256_set1_epi32(0);
    __m256i one = _mm256_set1_epi32(1);
    __m256d v = x;
    while (_mm256_testz_si256(n, zero))
    {
        x = _mm256_mul_pd(x, x);
        if (_mm256_testz_si256(n, one)) v = _mm256_mul_pd(v, x);
        ln = _mm_srli_si128(_mm256_extractf128_si256(n, 0), 1);
        hn = _mm_srli_si128(_mm256_extractf128_si256(n, 1), 1);
        n = _mm256_insertf128_si256(n, ln, 0);
        n = _mm256_insertf128_si256(n, hn, 1);
    }
    return v;
}
#endif


/* 
   _mm_srli_si128
   _mm256_insertf128_si256

   FICOM
   PCMPEQB
*/
