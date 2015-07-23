#ifndef _COMMON_TOOLS_H
#define _COMMON_TOOLS_H

#include <stdio.h>

FILE *open_if_present(const char *file, const char *mode);
int close_if_present(FILE *fh);

int fastq_type_to_offset(const char *type);

/* parse a comma-separated list of values into an array. check
   that the values are between 0 and 1 and increasing. */
void parse_csv_line(const char *csv, double *vals, unsigned *n_vals, unsigned max_vals);


#endif /* _COMMON_TOOLS_H */
