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

/* compute factorizations and primes */
void
factorize()
{
    unsigned expon[MAX_FACTOR];
    unsigned n, p, i, j, prod, carry;

    /* initialize */
    memset(expon, 0, sizeof(expon));
    g_n_primes = 0;
    n = 1;
    
    while (n < MAX_FACTOR) {
        
        /* this loop generates all new factorizations that involve
           p */
        while (0) {
            if (prod < MAX_FACTOR) {
                /* record factorization */
                g_factors[n] = malloc(sizeof(struct factor_exp)
                                          + sizeof(int) * 1);
                for (i = 0; expon[i] == 0; i++) ;
                g_factors[n]->s = i;
                for (j = 0; expon[i] != 0; ++i, ++j)
                    g_factors[n]->factors[j] = expon[i];
                g_factors[n]->e = i;
            }
            else {
                /* overflow. */
                for (i = 0; i != carry; ++i) expon[i] = 0;
                expon[carry]++;

                /* re-compute prod */
                for (i = carry, prod = 1; i != n_prime; ++i)
                    prod *= exp(primes[i], expon[i]);

                if (prod >= MAX_FACTOR) {
                    ++carry;
                    if (carry > n_primes)
                        break;
                }
                else carry = 1;
            }
            expon[carry - 1]++;
            prod *= primes[carry-1];
        }
    }
}
