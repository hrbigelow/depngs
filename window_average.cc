/*
  Input is a file containing <chromosome>, <position>, <depth> fields
  Options: -w <window_size> -e <every>
  Outputs to stdout a window-averaged depth, taking, for the current row,
  the window defined by p in (p - w/2, p + w/2)
  
  
 */

#include <getopt.h>
#include <cstdlib>
#include <cstdio>
#include <cstring>
// #include <cmath>
#include <vector>

#include "samutil/file_utils.h"
#include "histo.h"

struct ll
{
    int64_t pos;
    size_t depth;
    ll * next;
    ll(size_t _p, size_t _d, ll * _n) : pos(_p), depth(_d), next(_n) { }
};


int main(int argc, char *argv[])
{
    size_t window_size = 0, every = 1000, max_mem = 4e9, bins_per_unit = 1, nunits = 100;

    const char *label = "";
    const char *window_averaged_outfile = "";
    char c;
    while ((c = getopt(argc, argv, "w:e:l:m:b:n:a:")) >= 0)
    {
        switch(c)
        {
        case 'w': window_size = static_cast<size_t>(atof(optarg)); break;
        case 'e': every = static_cast<size_t>(atof(optarg)); break;
        case 'l': label = optarg; break;
        case 'm': max_mem = static_cast<size_t>(atof(optarg)); break;
        case 'b': bins_per_unit = static_cast<size_t>(atof(optarg)); break;
        case 'n': nunits = static_cast<size_t>(atof(optarg)); break;
        case 'a': window_averaged_outfile = optarg; break;
        }
    }

    size_t nbins = bins_per_unit * nunits;

    char const* infile = argv[optind];
    FILE * in_fh = fopen(infile, "r");
    if (in_fh == NULL)
    {
        fprintf(stderr, "Error: Couldn't open input file %s\n", infile);
        exit(10);
    }
    char const* histfile = argv[optind + 1];
    FILE *hist_fh = fopen(histfile, "w");
    if (! hist_fh)
    {
        fprintf(stderr, "Error: Couldn't histogram output file %s\n", histfile);
        exit(10);
    }
    FILE *win_fh = fopen(window_averaged_outfile, "w");
    
    char chrom[64], prev_chrom[64] = "";

    size_t depth, pileup_depth;
    int64_t pileup_pos = -1, window_pos = 1, half_window_size = window_size / 2;
    float running_sum = 0.0, current_count = 0.0;
    ll *head = NULL, *tail = NULL, *elem = NULL, *mid = NULL;

    size_t max_pileup_line_size = 1e7;
    size_t bytes_read, bytes_wanted = max_mem;

    char *chunk = new char[bytes_wanted + max_pileup_line_size + 1];
    size_t fread_nsec;
    char *last_fragment;

    size_t *hist = (size_t *)calloc(nbins, sizeof(size_t));
    size_t dbin;
    bool new_contig = true;

    //size_t line = 0;
    while (! feof(in_fh))
    {
        bytes_read = FileUtils::read_until_newline(chunk, bytes_wanted, max_pileup_line_size, in_fh, &fread_nsec);
        std::vector<char *> lines_vec = FileUtils::find_complete_lines_nullify(chunk, &last_fragment);
        size_t nlines = lines_vec.size();
        char **lines = lines_vec.data();
        char **line;

        for (line = lines; line != lines + nlines; )
        {
            if (pileup_pos < window_pos)
            {
                // we've gone past the pileup, need to load a new line.
                // this should happen whenever 
                sscanf(*line++, "%s\t%zi\t%*c\t%zu", chrom, &pileup_pos, &pileup_depth);
                new_contig = strcmp(chrom, prev_chrom);
            }
            else { new_contig = false; }
            if (new_contig) { window_pos = 1; }

            depth = pileup_pos == window_pos ? pileup_depth : 0;

            int64_t window_min_pos = window_pos - window_size;

            if (new_contig)
            {
                // starting a new contig.  flush remaining windows
                while (head != NULL && head->next != NULL)
                {
                    running_sum -= head->depth;
                    current_count--;
                    elem = head;
                    head = head->next;
                    while (mid != tail && mid->pos < head->pos + half_window_size)
                    {
                        mid = mid->next;
                    }
                    if ((tail->pos - head->pos) >= half_window_size)
                    {
                        if (win_fh && (mid->pos % every == 0))
                        {
                            fprintf(win_fh, "%s\t%s\t%zu\t%f\n", 
                                    label,
                                    prev_chrom, mid->pos, 
                                    running_sum / window_size);
                        }
                        dbin = bin(running_sum / window_size, bins_per_unit);

                        // tally histogram
                        if (dbin < nbins)
                        {
                            hist[dbin]++;
                        }

                    }
                    

                    delete elem;
                }

                if (strcmp(prev_chrom, "") != 0)
                {
                    print_hist(hist_fh, label, prev_chrom, hist, bins_per_unit, nbins);
                }

                // now start the new window
                head = new ll(window_pos, depth, NULL);
                tail = head;
                mid = head;
                current_count = 0;
                running_sum = 0;
            }
            else
            {
                // moving the window as usual.  append
                tail->next = new ll(window_pos, depth, NULL);
                running_sum += tail->depth;
                current_count++;
                tail = tail->next;

                // update mid until it is is within appropriate range of tail
                while (mid->pos < tail->pos - half_window_size)
                {
                    mid = mid->next;
                }

                // update head until it is within appropriate range of tail
                while (head->pos < window_min_pos)
                {
                    running_sum -= head->depth;
                    current_count--;
                    elem = head;
                    head = head->next;
                    delete elem;
                }
            }

            // output windows, update histogram
            if ((tail->pos - head->pos) >= half_window_size)
            {
                if (win_fh && (mid->pos % every == 0))
                {
                    fprintf(win_fh, "%s\t%s\t%zu\t%f\n", 
                            label,
                            prev_chrom, mid->pos, 
                            running_sum / window_size);
                }
                dbin = bin(running_sum / window_size, bins_per_unit);
                if (dbin < nbins)
                {
                    hist[dbin]++;
                }
            }
            strcpy(prev_chrom, chrom);
            ++window_pos;
        }
    }
    fclose(in_fh);

    // now we still have a non-empty linked list.  Proceed to print
    // out the remaining shrinking windows
    while (head->next != NULL)
    {
        running_sum -= head->depth;
        current_count--;
        elem = head;
        head = head->next;
        while (mid != tail && mid->pos < head->pos + half_window_size)
        {
            mid = mid->next;
        }
        if ((tail->pos - head->pos) >= half_window_size)
        {
            if (win_fh && (mid->pos % every == 0))
            {
                fprintf(win_fh, "%s\t%s\t%zu\t%f\n", 
                        label,
                        prev_chrom, mid->pos, 
                        running_sum / window_size);
            }
            dbin = bin(running_sum / window_size, bins_per_unit);
            if (dbin < nbins)
            {
                hist[dbin]++;
            }
        }
        delete elem;
    }
    print_hist(hist_fh, label, prev_chrom, hist, bins_per_unit, nbins);

    free(hist);
    fclose(hist_fh);
    if (win_fh) { fclose(win_fh); };
}
