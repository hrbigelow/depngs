/* functions for accurately computing mcn(a1,b1,c1) * mcn(a2,b2,c2) /
   mcn(a1+a2,b1+b2,c1+c2).  this is the multinomial adjustment factor
   for convolutions */

#define MAX_FACTOR 1000

/* will be populated with the primes, but only up to value of
   MAX_FACTOR, which will be less than capacity. */
static unsigned g_primes[MAX_FACTOR];

static unsigned g_n_primes;

/* represents a factorization of an integer in terms of the count of
   each prime number. factors holds the number of prime factors in
   [primes[s],...,primes[e-1]]*/
struct factor_exp {
    unsigned s, e;
    int factors[FLEX_ARRAY];
};

/* hold the factorizations for each number */
static struct factor_exp *g_factors[MAX_FACTOR];

/* factorial factorizations */
static struct factor_exp *g_factorial_pff[MAX_FACTOR];

/* compute factorizations and primes */
void
factorize()
{
    unsigned expon[MAX_FACTOR];
    unsigned n, p, s, e, i, j, prod, carry;

    /* initialize */
    memset(expon, 0, sizeof(expon));
    for (i = 0; i != MAX_FACTOR; ++i) g_factors[i] = NULL;

    g_n_primes = 0;
    n = 2;
    prod = 2;

    while (n != MAX_FACTOR) {
        primes[g_n_primes++] = n;

        /* this loop generates all new factorizations that involve
           p */
        while (0) {
            if (prod < MAX_FACTOR) {
                /* record factorization */
                for (s = 0; expon[s] == 0; s++) ;
                for (e = s + 1; expon[e] != 0; ++e) ;
                struct factor_exp *spp = malloc(sizeof(struct factor_exp) + sizeof(int) * (e - s));

                spp->s = s;
                spp->e = e;
                for (i = s, j = 0; i != e; ++i, ++j)
                    spp->factors[j] = expon[i];
                g_factors[prod] = spp;

                ++expon[0];
                prod *= primes[0];
            }
            else {
                /* overflow. */
                i = 0;
                while (i != g_n_primes && prod >= MAX_FACTOR) {
                    prod /= exp(primes[i], expon[i]);
                    expon[i++] = 0;
                    ++expon[i];
                    prod *= primes[i];
                }
                if (i == g_n_primes) break;
            }
        }
        /* advance n to next prime < MAX_FACTOR if it exists */
        while (g_factors[n++] && n != MAX_FACTOR)
            ;
    }
}


void
factorial_factorize()
{
    unsigned n, s = MAX_FACTOR, e = 0;
    int rsum[MAX_FACTOR];
    memset(rsum, 0, sizeof(rsum));

    for (n = 2; n != MAX_FACTOR; ++n) {
        s = MIN(s, g_factors[n]->s);
        e = MAX(e, g_factors[n]->e);
        for (i = s; i != e; ++i)
            rsum[i] += g_factors[n]->factors[i];
        struct factor_exp *fe = malloc(sizeof(struct factor_exp) + sizeof(int) * (e - s));
        fe->s = s;
        fe->e = e;
        memcpy(fe->factors, rsum + s, (e - s) * sizeof(rsum[0]));
        g_factorial_pff[n] = fe;
    }
}


/* compute (a+b+c)!/a!b!c!, returning a representation in pf form.
   assume a >= b >= c. */
struct factor_exp *
d_choose_abc(unsigned a, unsigned b, unsigned c)
{
    int rsum[MAX_FACTOR];
    unsigned d = a + b + c;
    assert(d < MAX_FACTOR);
    struct factor_exp *pf = g_factorial_pff[d];
    unsigned s = pf->s, e = pf->e;
    memcpy(rsum + s, pf->factors, (e - s) * sizeof(int));

    unsigned i, j;
    pf = g_factorial_pff[a];
    for (i = pf->s, j = 0; i != pf->e; ++i, ++j) rsum[i] -= pf->factors[j];

    pf = g_factorial_pff[b];
    for (i = pf->s, j = 0; i != pf->e; ++i, ++j) rsum[i] -= pf->factors[j];

    pf = g_factorial_pff[c];
    for (i = pf->s, j = 0; i != pf->e; ++i, ++j) rsum[i] -= pf->factors[j];

    while (rsum[s] == 0) ++s;
    while (rsum[e-1] == 0) --e;

    pf = malloc(sizeof(struct factor_exp) + (e - s) * sizeof(int));
    pf->s = s;
    pf->e = e;
    i = s;
    j = 0;
    while (i != e) pf->factors[j++] = rsum[i++];

    return pf;
}




/* reify a pff.  */
double
reify_pff(int *exp, unsigned s, unsigned e)
{
    double val = 1.0;
    unsigned i;
    for (i = 0; i != pff->e - pff->s; ++i)
        val *= exp();
}


/* compute the full adjustment coefficient */
double
adjust_coeff(struct factor_exp *n1,
             struct factor_exp *n2,
             struct factor_exp *d)
{
    int rsum[MAX_FACTOR];
    memset(rsum, 0, sizeof(rsum));
    memcpy(rsum + n1->s, n1->factors, (n1->e - n1->s) * sizeof(int));
    unsigned i, j;
    for (i = n2->s, j = 0; i != n2->e; ++i, ++j) rsum[i] += n2->factors[j];
    for (i = d->s, j = 0; i != d->e; ++i, ++j) rsum[i] -= d->factors[j];

    return reify_pff(rsum + s, s, e);
}
