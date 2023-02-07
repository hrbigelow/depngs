# dep - Diversity Estimation Probabilistically

A tool to estimate nucleotide base composition in non-clonal samples from next-gen
sequence alignment

## Synopsis

    dep dist [options] samples.rdb sample_pairings.rdb ref.fasta
    dep pileup [options] samples.rdb ref.fasta out.pileup

# Introduction

I first implemented while at the Broad Institute and further developed at Amgen.  The
problem setting is as follows.  The replicating virus or cancer genome can mutate at
each division, with the result that a population is non-clonal.  Each position in the
genome potentially has a different selective or mutational pressure.  Taking the
population as a whole, we want to measure the fraction of the population having A, C,
G, and T at each position in the genome.  Knowing this allows biologists to measure
the effect of a drug on the mutational profile, or to spot locations that correlate
with survival or death.   

While cancer cells divide and mutate, the immune system recognizes some of them
as bad, and is able to kill them off.  Other cells acquire mutations that allow
them to evade attack for awhile until the immune system adapts.  A Darwinian
natural selection dynamic takes hold.  In between adaptations, a mostly
decimated population of cancer cells can regrow from a very small subpopulation
('clonal expansion'). 

Think of the game like this:  The immune system is AlphaGo playing a billion Go
games simultaneously, each against a different amateur.  At each time step, the
amateurs each make a random move.  AlphaGo responds and the process repeats.
Also, at each time step, some fraction of games end, with the amateurs losing.
When these games end, there is room for any remaining amateurs to clone
themselves with the same board state until there are a billion games in play. 

So, at any given time, in a very rare cases, the amateur might make a really
good move, and since he can clone his board state, the overall setup provides a
brute-force opportunity for the amateurs to win.

Breaking from the metaphor now, the goal is to spot which positions in the
genome have undergone clonal expansion, and when.  This is done by comparing
the measured base composition of different cancer tissue samples at successive
timepoints.  And to get an accurate picture, one must be able to detect rare
subpopulations.

So the problem is this:  Now imagine I give you an urn with a billion balls,
some unknown fraction of red, blue, green, and yellow.  You can only take out
5000 balls, and you actually can't observe the color of each ball directly.
You have a machine that can measure the color, but it is not 100% accurate.  In
fact, it outputs a color and a confidence score (based on its own internal
metrics) telling you how likely it is to be correct.  For the sake of argument,
let's say we can take this confidence score as accurate.

So you get your 5000 measurements of (color, confidence score), and now you
want to estimate the fraction of different colors in the full urn.  The model I chose
is a two-stage model in which: 

# Implementation

dep is a multi-threaded filter which consumes one or more .bam files and outputs one
or more result files.  It is a unix-like filter in that it outputs an analysis for
each genomic position in the order given in the inputs.  The analysis can be any of:

- base composition estimation with confidence intervals
- mutational distance estimation with confidence intervals
- pileup

At the top level of execution it is a parallel for-loop with fixed memory buffer
specified by the user.  writing the output in the same order as the input, under the
race condition of different worker finishing times, presents a technical challenge.
This is solved using a linked list of output result buffers.  The linked-list
structure preserves the order.  Each buffer can be loaded by a worker thread, or
unloaded by a single writer thread.  To prevent worker thread starvation, there are
more buffers than threads.


