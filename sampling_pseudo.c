/*
A 'quantile' is a chosen mass partition fraction of a density.  A
'quantile value' is the position in the domain of the density where
that mass partition fraction is obtained.

Motivation:

Our distributions fall into 3 groups based on a quantile value qv of
interest.  qv is typically around 0.2 and represents the minimum
amount of non-reference base needed to make a locus 'interesting'.

Q is the 'quantile' function which provides the quantile corresponding
to a given quantile value of a distribution.  Otherwise known as the
CDF or Error Function.

A) non-changed             0.99 < Q(qv;I)
B) too-little-data         0.01 < Q(qv;I) < 0.99.
C) changed                 Q(qv;I) < 0.01

However, there are a few layers of indirection here.  First, we are
interested in P(*), not I(*).  Second, we can only get lower and
upper-bound estimates of Q(qv;I(*)), and these estimates converge
towards each other as we take more and more sample points.

We thus have:

bqv_lo, the 99% confidence lower bound estimate of Q(qv;I)
bqv_hi, the 99% confidence upper bound esitmate of Q(qv;I)

bqv_lo_adj, the 99% confidence lower bound esimtate of Q(qv;P)
bqv_hi_adj, the 99% confidence upper bound esimtate of Q(qv;P)

This then yields these three cases.  The vertical bars represent our
desired confidence in claiming a locus has 'less than qv non-reference
base' or 'greater than qv non-reference base'.  The brackets represent
our confidence in the lower and upper-bound estimates of Q(qv;P).

We have six cases:

 changed  | too-little-data  |  non-changed
I       [ |                  | ]
A         |                [ | ]
B       [ | ]                |
C         |                  |   [   ]
D         |       [   ]      |
E  [   ]  |                  |

I: More sample points are needed to say anything.  This is the most
   likely starting point when very few sample points have been taken.
   It should never happen that we end up here after taking max points.

A: 'changed' category eliminated.  More sample points needed to call
   'non-changed' or 'too-little-data'. If max points taken, do not
   report.  stop.

B: 'non-changed' category eliminated.  More sample points needed to
   call 'changed' or 'too-little-data'.  If max points taken, do not
   report. stop.

C: called as 'non-changed'.  stop.  (most common case)

D: called as 'too-little-data'.  stop.  (second most common case)

E: called as 'changed'.  Continue taking samples until max in order to
   get accurate estimates.  (Should there be a separate max for this?)

During a given loop, the possible paths through these states would be:

IBC or IC   (non-changed locus)
IAD or ID   (low-coverage locus that has no non-reference bases)
IBE or IE   (changed locus)
IB          (borderline non-changed locus, hitting max points)
IC          (borderline changed locus, hitting max points)


CONCEPTUAL NOTES:

We keep a tally of the number of loci that fall into each of the 6
categories.

Conceptually, our P distribution is wholly determined by the existing
data and thus can only represent changed, too-little-data, or
non-changed.  These three categories are determined solely by the
Q(qv;P) value being below lo_thresh, between lo_thresh and hi_thresh,
or above hi_thresh, respectively.  However, since we can't compute
Q(qv;P) directly, (can only generate confidence intervals for it)
cases D and E exist.  Faced with this, we choose to resolve the
sampling ambiguity in favor of the too-little-data category.

Loci whose true Q(qv;P) are in the changed partition but happen to
fall in case D will not be reported, but the number of such loci will
be tallied.  These loci represent the true borderline cases of
interest.  However, they are not reported, because one can always
adjust the threshold a little lower and get them reported.

In summary, no matter what choices of value or confidence theshold a
user chooses, there will be edge cases that question the threshold.
There is no way around it.

QV is the quantile value function.  It accepts a quantile and outputs
a quantile value, relative to a distribution function.
*/

/* 1. bqv is a 99% confidence lower bound estimate of Q(qv,I) */
float bqv_lo = QV(0.01, Beta(s + 1/2, n - s + 1/2));

if (s == 0 || s == n)
 {
     /* no useful information for applying weights */
     bqv_lo_adj = bqv;
 }
 else
 {
     w_hi = sum_i{P(x_i)/I(x_i)} s.t. x_i > qv;
     w_lo = sum_i{P(x_i)/I(x_i)} s.t. x_i < qv;
     w = w_hi + w_lo;
     c_hi = w_hi / s;
     c_lo = w_lo / (n - s);
     bqv_lo_adj = (c_hi * bqv) / (c_lo * (1 - bqv) + (c_hi * bqv));
 }

if (bqv_lo_adj > hi_thresh) break;
/* stop.  this is an obvious non-change. */
if (bqv_hi_adj < lo_thresh)


 else
   {
       /* generate another sample point and repeat from step 2. */
       /* once we hit the limiting number of points (1000), switch to
          'optimistic' phase. */      
   }

/* 5. The optimistic phase. */


