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
void parse_csv_line(const char *csv, double *vals, unsigned *n_vals)
{
    *n_vals = 0;
    errno = 0;
    double pv = -1;
    unsigned err = 0;
    char *loc;
    while (*csv != '\0' && *n_vals != MAX_NUM_QUANTILES)
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
                __func__, MAX_NUM_QUANTILES, csv);
        exit(1);
    }
}
