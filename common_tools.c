#include "common_tools.h"

#include <string.h>
#include <stdlib.h>
#include <errno.h>

// if file is not NULL and not '/dev/null', attempts to
// open the file.  In this case, it is an error if it cannot open the file.
// otherwise, returns NULL
FILE *open_if_present(const char *file, const char *mode)
{
    if (file == NULL 
        || strcmp(file, "/dev/null") == 0
        || strcmp(file, "") == 0)
        return NULL;
    else
    {
        FILE * fh = fopen(file, mode);
        if (fh == NULL)
        {
            fprintf(stderr, "Error: open_if_present: file %s not blank or '/dev/null' but"
                    " still couldn't open it\n", file);
            exit(1);
        }
        else return fh;
    }
}

int close_if_present(FILE *fh)
{
    if (fh != NULL) return fclose(fh);
    else return 0;
}


int fastq_type_to_offset(const char *type)
{
    static struct { const char *type; int offset; } 
    fastq_offsets[] = { 
        { "Sanger", 33 },
        { "Illumina18", 33 },
        { "Solexa", 64 },
        { "Illumina13", 64 },
        { "Illumina15", 64 }
    };

    unsigned i;
    for (i = 0; i != sizeof(fastq_offsets) / sizeof(fastq_offsets[0]); ++i)
        if (! strcmp(type, fastq_offsets[i].type)) return fastq_offsets[i].offset;
    return -1;
}


/* parse a comma-separated list of values into an array. check
   that the values are between 0 and 1 and increasing. */
void
parse_csv_line(const char *csv, double *vals, unsigned *n_vals, unsigned max_vals)
{
    *n_vals = 0;
    errno = 0;
    double pv = -1;
    unsigned err = 0;
    char *loc;
    while (*csv != '\0' && *n_vals != max_vals)
    {    
        if (*csv == ',') ++csv;
        *vals = strtod(csv, &loc);
        csv = loc;
        err = (err || *vals < 0 || *vals > 1 || *vals < pv);
        pv = *vals;
        ++vals;
        ++(*n_vals);
    }
    err = (err || *csv != '\0' || *n_vals == 0);
    if (errno || err)
    {
        fprintf(stderr, "%s: The input is supposed to be a comma-separated"
                "list of between 1 and %u increasing numbers in the [0, 1] range\n"
                "Input was: %s\n",
                __func__, max_vals, csv);
        exit(1);
    }
}


/* attempts to convert nptr to a double value using strtod.  if strtod
   doesn't consume exactly the entire contents of nptr, prints a
   message to standard error using conv_name as a tag, and then
   exits. */
double strtod_errmsg(const char *nptr, const char *conv_name)
{
    char *end;
    double rv = strtod(nptr, &end);
    if (*nptr == '\0' || end != nptr + strlen(nptr))
    {
        fprintf(stderr, "Error: %s: couldn't convert %s with value \"%s\" to a double\n",
                __func__, conv_name, nptr);
        exit(1);
    }
    return rv;
}


long int strtol_errmsg(const char *nptr, const char *conv_name)
{
    char *end;
    long int rv = strtol(nptr, &end, 10);
    if (*nptr == '\0' || end != nptr + strlen(nptr))
    {
        fprintf(stderr, "Error: %s: couldn't convert %s with value \"%s\" to a long int\n",
                __func__, conv_name, nptr);
        exit(1);
    }
    return rv;
}
