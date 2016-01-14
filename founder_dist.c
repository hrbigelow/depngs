/* A set of functions for calculating the distribution of the hidden
   founder base distribution from a set of basecalls and binned
   quality scores. */

/* compute the likelihood of hidden state composition with 'hs'
   secondary symbols, where the obsered configuration has op primary
   symbols and os secondary symbols, and the 'error' (miscall)
   probability is e. */
float
path_likelihood(unsigned op, unsigned os, unsigned hs, float e)
{
    float marg = 0;
    int p = hs > op ? op : hs;
    int s = hs > op ? hs - op : 0;
    for (; p >= 0 && s <= os; --p, ++s) {
        marg += m_choose_n(op, p)
            * m_choose_n(os, s)
            * exp(e, op + os - p - s)
            * exp(1 - e, p + s);
    }
    return marg;
}


/* find the mode of a component path function F(x;op;e) =
   m_choose_n(op,x) * e^(op-x)(1-e)^x using the ratio technique. */
unsigned
component_mode(unsigned op, float e)
{
    /* technique: Note that F(x+1)/F(x) = {(x+1)/(op-x)} * {(1-e)/e}.
       Thus, find the value of x such that (x+1)/(op-x) falls below
       (1-e)/e. */
    float fac = (1.0-e)/e;
    unsigned xp1, dif; /* 'x+1', 'op-x' */
    for (xp1 = 1, dif = op; xp1 != op; ++xp1, --dif) {
        if (xp1 / dif < fac) break;
    }
    return xp1 - 1;
}


/* compute the majority set of the conditional distribution over the
   composition of H, Cond(H;op,os,e) */
void
hidden_cond_dist(unsigned op, unsigned os, float e,
                 float **dist, unsigned *n_dist)
{
    /* find the component mode */
    unsigned hs = component_mode(op, e);
    unsigned end = hs + os;
    int s;
    unsigned p, s;
    float tot = 0;
    for (s = hs; s != end; ++s) {
        dist_tmp[s] = path_likelihood(op, os, s, e);
        tot += dist_tmp[s];
    }
    float thresh = tot * 1e-4; /* !!! parameter */
    float comp;
    /* search to the left and right */
    for (s = hs - 1; s >= 0, --s) {
        if ((comp = path_likelihood(op, os, s, e)) < thresh) break;
        dist_tmp[s] = comp;
    }
    for (s = end; s != op + os; ++s) {
        if ((comp = path_likelihood(op, os, s, e)) < thresh) break;
        dist_tmp[s] = comp;
    }
    
}

