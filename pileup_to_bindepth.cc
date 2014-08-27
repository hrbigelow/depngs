// Convert a pileup file to a bindepth format


#include <vector>
#include <string.h>

#include "samutil/file_utils.h"
#include "bindepth.h"

#define MIN(a,b) ((a) < (b) ? (a) : (b))


int main(int argc, char *argv[])
{
    if (argc != 5)
    {
        printf("\n%s <bufsize> <contig_dict_file> <pileup_input_file> <output_file>\n\n", argv[0]);
        printf("<bufsize> should be >= 1e7 for efficiently reading large pileup files\n\n");
        printf("<contig_dict_file> has the format:\n\n"
               "chr1\\t249250621\n"
               "chr2\\t243199373\n"
               "chr3\\t198022430\n"
               "...\n\n");
        exit(1);
    }
    size_t bufsize = static_cast<size_t>(atof(argv[1]));
    FILE *contig_dict_fh = fopen(argv[2], "r");
    FILE *in_fh = fopen(argv[3], "r");
    FILE *out_fh = fopen(argv[4], "w");

    if (! contig_dict_fh)
    {
        fprintf(stderr, "Error: couldn't open input file %s\n", argv[2]);
        exit(1);
    }
    if (! in_fh)
    {
        fprintf(stderr, "Error: couldn't open input file %s\n", argv[3]);
        exit(1);
    }
    if (! out_fh)
    {
        fprintf(stderr, "Error: couldn't open output file %s\n", argv[4]);
        exit(1);
    }

    // parse contig dictionary
    size_t ncontigs;
    contig_dict_t *contigs = parse_contig_dict(contig_dict_fh, &ncontigs);
    contig_dict_t *contig_def = contigs, *contig_def_end = contig_def + ncontigs;
    fclose(contig_dict_fh);

    // write contig dictionary to output file
    write_contig_dict(contigs, ncontigs, out_fh);

    size_t bytes_wanted = bufsize, bytes_read, max_pileup_line_size = 1e7, fread_nsec;
    char *chunk = (char *)malloc(sizeof(char) * (bufsize + max_pileup_line_size + 1));
    char *last_fragment;
    size_t pileup_pos = 0, locus_pos = 1, depth, pileup_depth;
    char contig[100], prev_contig[100] = "";
    bool new_contig = false;

    size_t npos = bufsize / sizeof(float);
    float *depth_buf = new float[npos];
    size_t locus_start_pos = 1;

    while (! feof(in_fh))
    {
        bytes_read = FileUtils::read_until_newline(chunk, bytes_wanted, max_pileup_line_size, in_fh, &fread_nsec);
        std::vector<char *> lines_vec = FileUtils::find_complete_lines_nullify(chunk, &last_fragment);
        size_t nlines = lines_vec.size();
        char **lines = lines_vec.data();
        char **line;

        for (line = lines; line != lines + nlines; )
        {
            if (pileup_pos < locus_pos)
            {
                sscanf(*line++, "%s\t%zu\t%*c\t%zu", contig, &pileup_pos, &pileup_depth);
                new_contig = strcmp(contig, prev_contig);
            }            
            else { new_contig = false; }

            if (new_contig)
            {
                if (strlen(prev_contig) != 0)
                {
                    // finish writing previous contig depth info
                    // locus_pos - 1 represents the last logical
                    // position on the previous contig.
                    // locus_pos - locus_start_pos represents the buffer index end
                    fwrite(depth_buf, sizeof(float), (locus_pos - locus_start_pos), out_fh);
                    size_t npos_remain = contig_def->size - (locus_pos - 1);

                    // zero the part of the buffer we're going to use
                    for (size_t f = 0; f != MIN(npos_remain, npos); ++f) { depth_buf[f] = 0l; }
                    while (npos_remain > 0)
                    {
                        fwrite(depth_buf, sizeof(float), MIN(npos_remain, npos), out_fh);
                        npos_remain -= MIN(npos_remain, npos);
                    }
                    locus_start_pos = 1;
                    
                    // alert the user of progress
                    fprintf(stderr, "Wrote %Zu loci depth for %s\n", contig_def->size, contig_def->name);
                    fflush(stderr);
                    ++contig_def;
                }
                

                // update contig_def further in case there are missing contigs
                while (strcmp(contig_def->name, contig) != 0 && contig_def != contig_def_end)
                {
                    // write zeros for any contigs that are completely missing
                    size_t npos_remain = contig_def->size;
                    for (size_t f = 0; f != MIN(npos_remain, npos); ++f) { depth_buf[f] = 0l; }
                    while (npos_remain > 0)
                    {
                        fwrite(depth_buf, sizeof(float), MIN(npos_remain, npos), out_fh);
                        npos_remain -= MIN(npos_remain, npos);
                    }
                    fprintf(stderr, "Wrote %Zu loci depth (all zeros) for %s\n", contig_def->size, contig_def->name);
                    fflush(stderr);

                    ++contig_def;
                }
                if (contig_def == contig_def_end)
                {
                    fprintf(stderr, "Error: Couldn't find contig %s in contig dictionary\n", contig);
                    exit(1);
                }

                // reset logical locus position
                locus_pos = 1;
            }

            depth = pileup_pos == locus_pos ? pileup_depth : 0;
            depth_buf[locus_pos++ - locus_start_pos] = depth;
            strcpy(prev_contig, contig);

            // signals that the buffer is full
            if (locus_pos - locus_start_pos == npos)
            {
                fwrite(depth_buf, sizeof(float), npos, out_fh);
                locus_start_pos = locus_pos;
            }
        }
    }

    // finish writing the remainder of the last contig.  locus_pos -
    //  locus_start_pos represents the buffer end.
    fwrite(depth_buf, sizeof(float), locus_pos - locus_start_pos, out_fh);

    size_t npos_remain = contig_def->size - (locus_pos - 1);
    for (size_t f = 0; f != MIN(npos_remain, npos); ++f) { depth_buf[f] = 0l; }
    while (npos_remain > 0)
    {
        fwrite(depth_buf, sizeof(float), MIN(npos_remain, npos), out_fh);
        npos_remain -= MIN(npos_remain, npos);
    }
    fprintf(stderr, "Wrote %Zu loci depth for %s\n", contig_def->size, contig_def->name);
    fflush(stderr);

    delete[] chunk;
    delete[] depth_buf;
    free(contigs);
    fclose(in_fh);
    fclose(out_fh);

    return 0;
}
