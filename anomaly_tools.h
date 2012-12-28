#ifndef _ANOMALY_TOOLS_H
#define _ANOMALY_TOOLS_H

class Posterior;
class NucleotideReader;
class LocusSummary;

#include "nucleotide_stats.h"

#include <cstddef>

double williams_moment_match_ratio(NucleotideStats const& stats,
                                   size_t N,
                                   double const* evaluation_point);


double strand_locus_anomaly_score(Posterior & posterior,
                                  JPD_DATA const& global_counts,
                                  LocusSummary const& full_locus,
                                  NucleotideReader const* data_reader,
                                  char strand,
                                  bool verbose);

double locus_anomaly_score(Posterior & posterior,
                           JPD_DATA const& global_counts,
                           LocusSummary const& full_locus,
                           NucleotideReader const* data_reader,
                                  bool verbose);


//compute relative entropy between P_global(b,q,s) and P_local(b,q,s)
//P_global(b,q,s) is computed from P_global(b,q,s,f) using a founder base (f)
//composition given as argument (typically the median)
double relative_entropy_anomaly(JPD_DATA const& global_counts,
                                LocusSummary const& locus_summary,
                                double est_composition[4]);

#endif // _ANOMALY_TOOLS_H
