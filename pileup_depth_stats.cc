/*
  Input a set of bindepth formatted pileup depth data, and a designation of
  which ones are to be used for calculating local averaged depth.

  Input a corresponding set of global sample-level average depth for each sample

  Input a window size, bins_per_unit, nunits.

  File
 */

#include "bindepth.h"
#include "histo.h"

#include <string.h>
#include <stdlib.h>
#include <getopt.h>
#include <sys/stat.h>
#include <assert.h>

char msg[] = 
    "Usage:\n\n"
    "pileup_depth_stats [options] <index_file> <output.hist>\n\n"
    "Options:\n\n"
    "-w  INT   window size for averaging depth [1]\n"
    "-m  INT   maximum memory to use in bytes [1e9]\n"
    "-b  INT   bins per unit-sized depth range for histogram binning [10]\n"
    "-n  INT   number of depth units for histogram [100]\n"
    "-a  STR   output file name for printing actual window-averaged depth at loci [empty]\n"
    "-e  INT   output averaged depth for every nth locus (relevant only with -a) [1]\n"
    "-l  STR   name of pseudo-sample representing local normalization [-1]\n"
    "-c  STR   name of a single selected contig to process.  If empty, process all contigs [empty]\n"
    "-r  STR   selected range of loci, i.e. 50000-10000.  If empty, process whole contigs [empty]\n"
    "-p  STR   parameter label string to attach to each output line. [na]\n"
    "-g  REAL  if local average is below this, locus is assigned 0 depth.  Should be in range [0,1].  [0.1]\n"
    "\n\n"
    "<index_file> has the format (tab-separated fields):\n\n"
    "<sample1_name>\t</path/to/sample1.bindepth>\t<use_as_normalizer(Y/N)>\t<global_average_depth>\t<space_delim_haploid_contigs>\n"
    "<sample2_name>\t</path/to/sample2.bindepth>\t<use_as_normalizer(Y/N)>\t<global_average_depth>\t<space_delim_haploid_contigs>\n"
    "...\n\n"
    "<use_as_normalizer(Y/N)> is used to mark those samples that are expected to have normal\n"
    "ploidy for all contigs.  These samples will be used for local normalization.\n\n"
    "<global_average_depth>: if this is the first run on these samples, this should be set to the value 1\n"
    "for every sample.  Then, subsequent runs can plug in the value at the peak of the histogram\n"
    "for each sample.\n\n"
    "<space_delim_haploid_contigs> is a space-delimited list of contig names that are expected\n"
    "to have haploid content. It is typically either 'none' signifying a female, or 'X Y' signifying\n"
    "a male sample (haploid X, haploid Y).\n\n";

#define MIN(a,b) ((a) < (b) ? (a) : (b))
#define MAX(a,b) ((a) > (b) ? (a) : (b))

struct file_and_avg_t
{
    FILE *fh;
    char name[20];
    bool normalizer;
    float avg_depth;
    char haploid_contigs[100];
};


size_t chunk_size, window_size = 1;


// since we do a memmmove, the last window contained in the previous chunk
// is also contained in the current chunk.  it is processed as the first window,
// unless this is the last chunk
void refresh_chunk(file_and_avg_t *sample,
                   size_t nsamples,
                   contig_dict_t *ctg,
                   int64_t spos,
                   int64_t epos,
                   float **d, float *wd,
                   int64_t *cpos, size_t *ce,
                   bool *last_chunk)
{
    size_t jump_size, read_size, overlap, region_size = (size_t)(epos - spos);
    if (*cpos == -1)
    {
        // indicates we are at the start of a new contig
        for (size_t s = 0; s != nsamples; ++s) { fseek(sample[s].fh, spos * sizeof(float), SEEK_CUR); }
        jump_size = 0;
        read_size = MIN(chunk_size, region_size);
        overlap = 0;
        *last_chunk = read_size == region_size;
        *cpos = spos;
        *ce = read_size - window_size + ((*last_chunk) ? 1 : 0);
    }
    else
    {
        jump_size = chunk_size - window_size;
        read_size = MIN(epos - (*cpos + *ce + window_size), jump_size);
        overlap = window_size;
        *last_chunk = (size_t)(epos - *cpos) <= chunk_size;
        *cpos += jump_size;
        *ce = read_size + ((*last_chunk) ? 1 : 0);
    }

    float *dp, *dpe, global_norm;
    for (size_t s = 0; s != nsamples; ++s) 
    { 
        // move the d by jump_size and initialize them.
        memmove(d[s], d[s] + jump_size, overlap * sizeof(float));
        fread(d[s] + overlap, sizeof(float), read_size, sample[s].fh);
        dp = d[s] + overlap, dpe = dp + read_size;
        global_norm = sample[s].avg_depth * window_size;
        while (dp != dpe) { *dp++ /= global_norm; }
        if (*last_chunk) { fseek(sample[s].fh, (ctg->size - epos) * sizeof(float), SEEK_CUR); }
    }

    // must compute the first window-averaged depth if this is the
    // first reading on this contig
    if (overlap == 0)
    {
        for (size_t s = 0; s != nsamples; ++s) 
        { 
            wd[s] = 0;
            for (dp = d[s]; dp != d[s] + window_size; ++dp) { wd[s] += *dp; }
        }
    }

}


file_and_avg_t *parse_index(FILE *in_fh, size_t *nsamples)
{
    struct stat s;
    size_t sz = (fstat(fileno(in_fh), &s), s.st_size);
    char *contents = (char *)malloc(sz), *c = contents, *cend = c + sz;
    fread(contents, sizeof(char), sz, in_fh);
    
    *nsamples = 0;
    while (c != cend) { ++(*nsamples); c = strchr(c, '\n') + 1; }

    file_and_avg_t *samples = (file_and_avg_t *)malloc(sizeof(file_and_avg_t) * *nsamples), *sp = samples;
    c = contents;
    char yn, path[512];
    while (c != cend)
    {
        sp->haploid_contigs[0] = ' ';
        sscanf(c, "%s\t%s\t%c\t%f\t%[^\t\n]\n", sp->name, path, &yn, &sp->avg_depth, sp->haploid_contigs + 1);
        strcat(sp->haploid_contigs, " ");
        sp->normalizer = (yn == 'Y');
        sp->fh = fopen(path, "r");
        c = strchr(c, '\n') + 1;
        ++sp;
    }
    free(contents);
    return samples;
}


int main(int argc, char *argv[])
{
    size_t every = 1000, max_mem = 4e9, bins_per_unit = 10, nunits = 100;
    int64_t selected_spos = -1, selected_epos = -1;
    float min_reliable_local_avg = 0.1;
    
    char *window_averaged_outfile = NULL;
    char *selected_contig = NULL, *selected_range = NULL;
    const char *local_norm_string = "-1";
    const char *param_label = "na";
    char c;
    while ((c = getopt(argc, argv, "w:e:m:b:n:a:l:c:r:p:g:")) >= 0)
    {
        switch(c)
        {
        case 'w': window_size = static_cast<size_t>(atof(optarg)); break;
        case 'e': every = static_cast<size_t>(atof(optarg)); break;
        case 'm': max_mem = static_cast<size_t>(atof(optarg)); break;
        case 'b': bins_per_unit = static_cast<size_t>(atof(optarg)); break;
        case 'n': nunits = static_cast<size_t>(atof(optarg)); break;
        case 'a': window_averaged_outfile = optarg; break;
        case 'l': local_norm_string = optarg; break;
        case 'c': selected_contig = optarg; break;
        case 'r': selected_range = optarg; break;
        case 'p': param_label = optarg; break;
        case 'g': min_reliable_local_avg = atof(optarg); break;
        default: fprintf(stderr, msg); return 1; break;
        }
    }

    if (argc != optind + 2)
    {
        fprintf(stderr, "%s", msg);
        return 1;
    }

    if (selected_range) { 
        int c;
        if ((c = sscanf(selected_range, "%zi-%zi", &selected_spos, &selected_epos)) != 2)
        {
            fprintf(stderr, "Error: selected range \"%s\" (argument -r) should be like '0-1000'\n", 
                    selected_range);
            exit(2);
        }
    }

    size_t nbins = bins_per_unit * nunits;
    size_t half_window_size = window_size / 2;

    const char *index_file  = argv[optind];
    FILE * in_fh = fopen(index_file, "r");
    if (in_fh == NULL)
    {
        fprintf(stderr, "Error: Couldn't open input file %s\n", index_file);
        exit(10);
    }
    
    size_t nsamples;
    file_and_avg_t *sample = parse_index(in_fh, &nsamples);

    char const* histfile = argv[optind + 1];
    FILE *hist_fh = fopen(histfile, "w");
    if (! hist_fh)
    {
        fprintf(stderr, "Error: Couldn't histogram output file %s\n", histfile);
        exit(10);
    }

    FILE *win_fh = fopen(window_averaged_outfile, "w");

    size_t nloci = max_mem / sizeof(float);
    chunk_size = nloci / nsamples;

    float *dbuf = new float[nloci], *dbufp;
    float **d = new float*[nsamples], **dp;

    size_t *hbuf = new size_t[nbins * (nsamples + 1)], *hbufp, *hbufe = hbuf + (nbins * (nsamples + 1));
    size_t **h = new size_t*[nsamples + 1], **hp;
    size_t s, ncontigs;

    // initialize files
    contig_dict_t *ctg, **contigs = new contig_dict_t*[nsamples];
    for (s = 0, dbufp = dbuf, hbufp = hbuf; s != nsamples; ++s)
    { 
        d[s] = dbufp; dbufp += chunk_size; 
        h[s] = hbufp; hbufp += nbins;
        contigs[s] = read_contig_dict(sample[s].fh, &ncontigs);
    }
    h[nsamples] = hbufp;

    // here, it would be prudent to check that all input files have
    // identical contig dictionaries.

    // main loop
    /*
      ==> At start of contig:
      0. cpos = -1 (cpos: position of first locus in chunk in contig)                                 
      1. read chunks d[N]                           
      2. d[N] /= (avg_depth * window_size)
      3. calcualate initial running averages wd[N]
      4. ci = 0 (ci: index into any chunk)

      ==> At each locus step, until the end of a chunk
      5. calculate A = avg(wd[N]) for 'normalizers'
      6. calculate wdn[N] as wd[N] / A
      7. tally histograms hist[N] from the wdn[N]
      8. if cpos + ci + hw % every == 0, output wdn[N]
      9. wd[N] += (d[ci + w] - d[ci])

      ==> At end of chunk: (ci == csize - w):
      1. memmove d[N][ci - w] values to the start
      2. read partial chunks of size csize - w
      3. d[N] /= (avg_depth * window_size)

     */

    float *wd = new float[nsamples], *wdp, *wde = wd + nsamples;
    float *wdn = new float[nsamples], *wdnp, *wdne = wdn + nsamples;
    float local_avg, **wd_norms = new float*[nsamples], **wd_normsp = wd_norms, **wd_normse;

    size_t nnormalizers = 0;
    for (s = 0; s != nsamples; ++s)
    {
        nnormalizers += sample[s].normalizer ? 1 : 0;
        if (sample[s].normalizer) { *wd_normsp = &wd[s]; ++wd_normsp; }
    }
    wd_normse = wd_norms + nnormalizers;

    // use the first contig dictionary as a representative
    size_t ci, ce, locus;
    int64_t cpos, spos, epos;
    bool last_chunk;
    
    float *ploidybuf = new float[nsamples];
    float *ploidybuf_p;
    char ploidy_query[100] = "";

    for (ctg = contigs[0]; ctg != contigs[0] + ncontigs; ++ctg)
    {
        cpos = -1; // indicates 'start of contig'
        local_avg = 0;

        ploidybuf_p = NULL;
        for (s = 0; s != nsamples; ++s)
        {
            strcpy(ploidy_query, " ");
            strcat(ploidy_query, ctg->shortname); // create a query that 
            strcat(ploidy_query, " ");

            ploidybuf[s] = strstr(sample[s].haploid_contigs, ploidy_query) 
                ? (ploidybuf_p = ploidybuf, 0.5) : 1.0;
        }

        // zero out the histogram buffer
        for (hbufp = hbuf; hbufp != hbufe; ++hbufp) { *hbufp = 0; }

        if (selected_contig && strcmp(ctg->name, selected_contig))
        {
            // skip this contig
            for (s = 0; s != nsamples; ++s)
            {
                fseek(sample[s].fh, ctg->size * sizeof(float), SEEK_CUR);
            }
            fprintf(stderr, "Skipped writing stats for %s (%Zu bases)\n", ctg->name, ctg->size);
            fflush(stderr);
            continue;
        }


        last_chunk = false;
        while (! last_chunk)
        {
            spos = selected_spos == -1 ? 0 : MIN((size_t)selected_spos, ctg->size);
            epos = selected_epos == -1 ? ctg->size : MIN((size_t)selected_epos, ctg->size);
            
            refresh_chunk(sample, nsamples, ctg, spos, epos, d, wd, &cpos, &ce, &last_chunk);

            // process all but the last window that fits in this chunk.
            // ce == chunk_size - window_size
            
            for (ci = 0; ci != ce; ++ci)
            {
                // compute local_avg normalizer.  if there are no normalizers, use the dummy value of 1
                if (wd_norms == wd_normse) { local_avg = 1; }
                else if (ploidybuf_p)
                {
                    // expected ploidy is mixed
                    local_avg = 0;
                    wd_normsp = wd_norms;
                    while (wd_normsp != wd_normse) { local_avg += **wd_normsp / *ploidybuf_p++; ++wd_normsp; }
                    local_avg /= nnormalizers;
                    ploidybuf_p = ploidybuf;
                }
                else
                {
                    local_avg = 0;
                    wd_normsp = wd_norms;
                    while (wd_normsp != wd_normse) { local_avg += **wd_normsp; ++wd_normsp; }
                    local_avg /= nnormalizers;
                }

                // calculate local averages, update histograms
                for (wdnp = wdn, wdp = wd, hp = h; wdnp != wdne; ++wdnp, ++wdp)
                {
                    // if the local_avg is unreliable, then only the extreme exceptions should
                    // be recorded
                    *wdnp = (local_avg < min_reliable_local_avg && *wdp < 0.5) ? 0 : MAX(*wdp / local_avg, 0);
                    (*hp++)[MIN(bin(*wdnp, bins_per_unit), nbins - 1)]++;
                }
                // calculate special histogram for local average
                (*hp)[MIN(bin(local_avg, bins_per_unit), nbins - 1)]++;

                // output selected depth lines
                if (win_fh && (locus = cpos + ci + half_window_size) % every == 0)
                {
                    for (s = 0; s != nsamples; ++s) 
                    { 
                        fprintf(win_fh, "%s\t%s\t%zu\t%6.4f\t%zu\t%s\n", 
                                sample[s].name, ctg->shortname, locus, wdn[s], window_size, param_label);
                    }
                    fprintf(win_fh, "%s\t%s\t%zu\t%6.4f\t%zu\t%s\n", 
                            local_norm_string, ctg->shortname, locus, local_avg, window_size, param_label);
                }

                // update wd, moving it forward by one
                for (wdp = wd, dp = d; wdp != wde; ++wdp, ++dp)
                { 
                    *wdp += ((*dp)[ci + window_size] - (*dp)[ci]);
                    // assert(*wdp < 1e20);
                }
            }
        }        

        // print out the histograms
        for (s = 0; s != nsamples; ++s) 
        {
            for (size_t bin = 0; bin != nbins; ++bin)
            {
                fprintf(hist_fh, "%s\t%s\t%6.4f\t%Zu\t%Zu\t%s\n",
                        sample[s].name, ctg->shortname, 
                        bin_center(bin, bins_per_unit), h[s][bin],
                        window_size, param_label);
            }
        }
        // output special histogram
        for (size_t bin = 0; bin != nbins; ++bin)
        {
            fprintf(hist_fh, "%s\t%s\t%6.4f\t%Zu\t%Zu\t%s\n",
                    local_norm_string, ctg->shortname, 
                    bin_center(bin, bins_per_unit), h[nsamples][bin],
                    window_size, param_label);
        }
        fflush(hist_fh);

        fprintf(stderr, "Finished stats for %s (%Zu bases)\n", ctg->name, ctg->size);
        fflush(stderr);
    }

    if (win_fh) { fclose(win_fh); }

    fclose(hist_fh);

    for (size_t f = 0; f != nsamples; ++f)
    {
        if (sample[f].fh) { fclose(sample[f].fh); }
        free(contigs[f]);
    }
    delete[] contigs;
    delete[] sample;
    delete[] dbuf;
    delete[] d;
    delete[] hbuf;
    delete[] h;
    delete[] wd;
    delete[] wdn;
    delete[] wd_norms;
    delete[] ploidybuf;
}
