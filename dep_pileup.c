/* Convert one BAM file to a pileup file, using a streaming,
   multi-threaded approach. */
#include "bam_sample_info.h"


int pileup_usage()
{
    char *tmp_require = bam_flag2str(mplp->rflag_require);
    char *tmp_filter  = bam_flag2str(mplp->rflag_filter);
    
    // Display usage information, formatted for the standard 80 columns.
    // (The unusual string formatting here aids the readability of this
    // source code in 80 columns, to the extent that's possible.)

    fprintf(stderr,
            "\nUsage: dep pileup [options] in.fasta in.bam [in2.bam ...]\n\n"
            "Options:\n"
            "-t  INT    number of threads [1]\n"
            "-F  INT    Output Phred Quality Encoding offset (33 or 64) [33]\n"
            "-A  FLAG   do not discard anomalous read pairs\n"
            "-B  FLAG   disable BAQ (per-Base Alignment Quality)\n"
            "-C  INT    adjust mapping quality; recommended:50, disable:0 [0]\n"
            "-E  FLAG   recalculate BAQ on the fly, ignore existing BQs\n"
            "-G  STR    file list of read groups to exclude\n"
            "-l  STR    file with ranges (chr start end) to process\n"
            "-q  INT    min mapQ to include alignments [%d]\n"
            "-Q  INT    min baseQ/BAQ to include in pileup [%d]\n"
            "-x  STR    flags that must be present for inclusion (string or int) [%s]\n"
            "-X  STR    flags that must be absent for inclusion (string or int) [%s]\n",
            mplp->min_mq,
            mplp->min_baseQ,
            tmp_require,
            tmp_filter);
    
    free(tmp_require);
    free(tmp_filter);
}

int main_pileup(int argc, char **argv)
{
    char c;
    
    /* adapted from samtools/bam_plcmd.c.*/
    while ((c = getopt(argc, argv, "t:F:ABC:EG:l:q:Q:x:X:")) >= 0) {
        switch (c) {
        case 't': n_threads = strtol_errmsg(optarg, "-t (n_threads)"); break;
        case 'F': phred_offset = strtol_errmsg(optarg, "-F (phred_offset)"); break;
        case 'A': use_orphan = 1; break;
        case 'B': mplp.flag &= ~MPLP_REALN; break;
        case 'C': mplp.capQ_thres = atoi(optarg); break;
        case 'E': mplp.flag |= MPLP_REDO_BAQ; break;
        case 'G': {
                FILE *fp_rg;
                char buf[1024];
                mplp.rghash = khash_str2int_init();
                if ((fp_rg = fopen(optarg, "r")) == 0)
                    fprintf(stderr, "(%s) Fail to open file %s. Continue anyway.\n", __func__, optarg);
                while (!feof(fp_rg) && fscanf(fp_rg, "%s", buf) > 0) // this is not a good style, but forgive me...
                    khash_str2int_inc(mplp.rghash, strdup(buf));
                fclose(fp_rg);
            }
            break;
        case 'l':
                  mplp.bed = bed_read(optarg);
                  if (!mplp.bed) { print_error_errno("Could not read file \"%s\"", optarg); return 1; }
                  break;
        case 'q': mplp.min_mq = atoi(optarg); break;
        case 'Q': mplp.min_baseQ = atoi(optarg); break;

        case 'x' :
            mplp.rflag_require = bam_str2flag(optarg);
            if ( mplp.rflag_require<0 ) { fprintf(stderr,"Could not parse -x %s\n", optarg); return 1; }
            break;
        case 'X' :
            mplp.rflag_filter = bam_str2flag(optarg);
            if ( mplp.rflag_filter<0 ) { fprintf(stderr,"Could not parse --X %s\n", optarg); return 1; }
            break;

        case  4 : mplp.openQ = atoi(optarg); break;
        case 'p': mplp.flag |= MPLP_PER_SAMPLE; break;
        case 'D': mplp.fmt_flag |= B2B_FMT_DP; fprintf(stderr, "[warning] samtools mpileup option `-D` is functional, but deprecated. Please switch to `-t DP` in future.\n"); break;
        case 'S': mplp.fmt_flag |= B2B_FMT_SP; fprintf(stderr, "[warning] samtools mpileup option `-S` is functional, but deprecated. Please switch to `-t SP` in future.\n"); break;
        case 'V': mplp.fmt_flag |= B2B_FMT_DV; fprintf(stderr, "[warning] samtools mpileup option `-V` is functional, but deprecated. Please switch to `-t DV` in future.\n"); break;

        default:
            fprintf(stderr,"Invalid option: '%c'\n", c);
            return 1;
        }
    }

    if (argc - optind != 3) return pileup_usage();
        
    char *fasta_file = argv[optind];
    char *bam_file = argv[optind + 1];
    char *plp_file = argv[optind + 2];

    struct thread_queue *tq = 
        pileup_init(samples_file,
                    fasta_file,
                    query_range_file,
                    n_threads,
                    n_readers,
                    max_input_mem,
                    min_quality);

    printf("Starting input processing.\n");
    thread_queue_run(tq);
    thread_queue_free(tq);

    printf("Finished.\n");
    
    return 0;
}


static struct {
    struct bam_reader_par *reader_buf;
    void **reader_par;
    unsigned n_readers;
    struct pair_ordering_range *ranges;
    unsigned n_ranges;
    unsigned n_threads;
} thread_params;



/* provide thread_queue API functions */
void
pileup_worker(const struct managed_buf *in_bufs,
              struct managed_buf *out_bufs)
{
    struct managed_buf bam = { NULL, 0, 0 };
    char **out = malloc(bam_samples.n * sizeof(char *));
    unsigned s;
    for (s = 0; s != bam_samples.n; ++s) {
        bam_inflate(&in_bufs[s], &bam);
        tally_pileup_stats(bam, s);
        summarize_pileup_stats(s);
        out[s] = out_bufs[s].buf;
        out_bufs[s].size = 0;
    }
    struct pileup_data pdat;
    struct pileup_locus_info ploc;
    unsigned add, more_loci = 1;
    while (more_loci) {
        pileup_current_info(&ploc);
        for (s = 0; s != bam_samples.n; ++s) {
            pileup_current_data(s, &pdat);
            add = pdat.calls.size + pdat.quals.size + MAX_LABEL_LEN + 5;
            out_bufs[s].size += add;
            ALLOC_GROW_REMAP(out_bufs[s].buf, out[s], 
                             out_bufs[s].size,
                             out_bufs[s].alloc);
            out[s] += 
                sprintf(out[s], "%s\t%s\t%u\t%c\t%u\t",
                        bam_samples.m[s].label,
                        ploc.refname,
                        ploc.pos,
                        ploc.refbase,
                        pdat.n_match_hi_q);
            
            strncpy(out[s], pdat.calls.buf, pdat.calls.size);
            out[s] += pdat.calls.size;
            *out[s]++ = '\t';
            strncpy(out[s], pdat.quals.buf, pdat.quals.size);
            *out[s]++ = '\n';
        }
        more_loci = pileup_next_pos();
    }   

    /* frees statistics that have already been used in one of the
       distance calculations. */
    pileup_clear_finished_stats();
    free(out);
}


/* just write to stdout.  there is only one output file */
void
pileup_offload(void *par, const struct managed_buf *bufs)
{
    write(0, bufs[0].buf, bufs[0].size);
}

void
pileup_oncreate()
{
    batch_pileup_thread_init(n_samples);
}


void
pileup_onexit()
{
    batch_pileup_thread_free();
}


struct thread_queue *
pileup_init(const char *samples_file,
            const char *fasta_file,
            const char *query_range_file,
            unsigned n_threads,
            unsigned n_readers,
            unsigned long max_input_mem,
            unsigned min_qual)
{
    bam_sample_info_init(samples_file, NULL);
    batch_pileup_init(min_qual, fasta_file);

    unsigned long n_total_loci;
    thread_params.ranges = 
        parse_query_ranges(query_range_file,
                           &thread_params.n_ranges,
                           &n_total_loci);

    thread_params.reader_par = (void **)malloc(n_readers * sizeof(void *));

    unsigned r, s;
    for (r = 0; r != n_readers; ++r) {
        thread_params.reader_buf[r] = (struct bam_reader_par){
            malloc(bam_samples.n * sizeof(struct bam_stats)),
            bam_samples.n,
            thread_params.ranges, 
            thread_params.ranges + thread_params.n_ranges
        };
        for (s = 0; s != bam_samples.n; ++s)
            bam_stats_init(bam_samples.m[s].bam_file, 
                           &thread_params.reader_buf[r].m[s]);
        
        thread_params.reader_par[r] = &thread_params.reader_buf[r];
    }

    /* chunking strategy */
    if (query_range_file)
        cs_init_by_range(n_total_loci, bam_samples.n);
    else
        cs_init_whole_file(bam_samples.n);

#define MAX_BYTES_SMALL_CHUNK 1e9
#define SMALL_CHUNK 5e6
#define DEFAULT_BYTES_PER_LOCUS 100

    cs_set_defaults(MAX_BYTES_SMALL_CHUNK,
                    SMALL_CHUNK, 
                    DEFAULT_BYTES_PER_LOCUS);

    unsigned n_extra = n_threads * 2;
    unsigned n_outputs = 1; /* just producing a pileup file */

    struct thread_queue_reader_t reader = { bam_reader, bam_scanner };
    struct thread_queue *tq =
        thread_queue_init(reader, reader_par,
                          pileup_worker,
                          pileup_offload, NULL,
                          pileup_on_create,
                          pileup_on_exit,
                          n_threads, n_extra, n_readers, bam_samples.n,
                          1, max_input_mem);

    return tq;
}


void
pileup_free()
{
    bam_sample_info_free();
    batch_pileup_free();

    unsigned r, s;
    for (r = 0; r != thread_params.n_readers; ++r) {
        for (s = 0; bam_samples.n; ++s)
            bam_stats_free(&thread_params.reader_buf[r].m[s]);
        free(thread_params.reader_buf[r].m);
    }

    free(thread_params.reader_buf);
    free(thread_params.reader_par);
    free(thread_params.ranges);
}
