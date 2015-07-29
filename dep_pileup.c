/* Convert one BAM file to a pileup file, using a streaming,
   multi-threaded approach. */
#include "bam_sample_info.h"
#include <stdlib.h>
#include "htslib/sam.h"
#include "khash.h"
#include "thread_queue.h"
#include "chunk_strategy.h"
#include "batch_pileup.h"
#include "bam_reader.h"
#include "common_tools.h"

static struct {
    unsigned min_base_quality; /* minimum base quality to be included in pileup */
    unsigned min_map_quality;
    unsigned long max_mem;
    int rflag_require;
    int rflag_filter;
    unsigned n_threads;
    unsigned n_readers;
    unsigned phred_offset;
} opts = { 5, 10, 1e9, 0, 0, 1, 1, 33 };


/* a place to hold allocated input for thread_queue that thread_queue
   does not own. */
static struct {
    struct bam_reader_par *reader_buf;
    void **reader_par;
    unsigned n_readers;
    struct pair_ordering_range *ranges;
    unsigned n_ranges;
    FILE *pileup_fh;
} thread_params;



int
pileup_usage()
{
    char *tmp_require = bam_flag2str(opts.rflag_require);
    char *tmp_filter  = bam_flag2str(opts.rflag_filter);

    fprintf(stderr,
            "\nUsage: dep pileup [options] samples.rdb in.fasta out.pileup\n\n"
            "Options:\n"
            "-t  INT    number of threads [%d]\n"
            "-R  INT    number of concurrent readers to allow (<= -t option) [%d]\n"
            "-m  INT    number bytes of memory to use [%ld]\n"
            "-F  INT    Output Phred Quality Encoding offset (33 or 64) [%d]\n"
            "-A  FLAG   do not discard anomalous read pairs\n"
            "-B  FLAG   disable BAQ (per-Base Alignment Quality)\n"
            "-C  INT    adjust mapping quality; recommended:50, disable:0 [0]\n"
            "-E  FLAG   recalculate BAQ on the fly, ignore existing BQs\n"
            "-G  STR    file list of read groups to exclude\n"
            "-l  STR    file with locus ranges (chr start end) to process [empty = process all input]\n"
            "-q  INT    min mapping quality to include alignments [%d]\n"
            "-Q  INT    min base_quality/BAQ to include in pileup [%d]\n"
            "-x  STR    flags that must be present for inclusion (string or int) [%s]\n"
            "-X  STR    flags that must be absent for inclusion (string or int) [%s]\n"
            "\n\n"
            "samples.rdb has lines of <sample_label>\t</path/to/sample.bam>\n\n"
            "There must also exist an <in.fasta>.fai file produced by samtools faidx\n\n",
            opts.n_threads,
            opts.n_readers,
            opts.max_mem,
            opts.phred_offset,
            opts.min_map_quality,
            opts.min_base_quality,
            tmp_require,
            tmp_filter);
    
    free(tmp_require);
    free(tmp_filter);
    return 1;
}


struct thread_queue *
pileup_init(const char *samples_file,
            const char *fasta_file,
            const char *pileup_file,
            const char *query_range_file,
            unsigned n_threads,
            unsigned n_readers,
            unsigned long max_input_mem,
            unsigned min_qual);

/* extern char *optarg; */
/* extern int optind; */


int main_pileup(int argc, char **argv)
{
    int c;
    char *query_range_file = NULL;

    /* adapted from samtools/bam_plcmd.c.*/
    while ((c = getopt(argc, argv, "t:R:m:F:ABC:EG:l:q:Q:x:X:")) >= 0) {
        switch (c) {
        case 't': opts.n_threads = strtol_errmsg(optarg, "-t (n_threads)"); break;
        case 'm': opts.max_mem = (size_t)strtod_errmsg(optarg, "-m (max_mem)"); break; 
        case 'R': opts.n_readers = strtol_errmsg(optarg, "-R (n_readers)"); break;
        case 'F': opts.phred_offset = strtol_errmsg(optarg, "-F (phred_offset)"); break;
        /* case 'A': opts.use_orphan = 1; break; */
        /* case 'B': mplp.flag &= ~MPLP_REALN; break; */
        /* case 'C': mplp.capQ_thres = atoi(optarg); break; */
        /* case 'E': mplp.flag |= MPLP_REDO_BAQ; break; */
        /* case 'G': { */
        /*         FILE *fp_rg; */
        /*         char buf[1024]; */
        /*         mplp.rghash = khash_str2int_init(); */
        /*         if ((fp_rg = fopen(optarg, "r")) == 0) */
        /*             fprintf(stderr, "(%s) Fail to open file %s. Continue anyway.\n", __func__, optarg); */
        /*         while (!feof(fp_rg) && fscanf(fp_rg, "%s", buf) > 0) // this is not a good style, but forgive me... */
        /*             khash_str2int_inc(mplp.rghash, strdup(buf)); */
        /*         fclose(fp_rg); */
        /*     } */
        /*     break; */
        case 'l': query_range_file = optarg; break;
        case 'q': opts.min_map_quality = strtol_errmsg(optarg, "-q (min_map_quality)"); break;
        case 'Q': opts.min_base_quality = strtol_errmsg(optarg, "-Q (min_base_quality)"); break;

        case 'x' :
            opts.rflag_require = bam_str2flag(optarg);
            if ( opts.rflag_require<0 ) { fprintf(stderr,"Could not parse -x %s\n", optarg); return 1; }
            break;
        case 'X' :
            opts.rflag_filter = bam_str2flag(optarg);
            if ( opts.rflag_filter<0 ) { fprintf(stderr,"Could not parse --X %s\n", optarg); return 1; }
            break;

        /* case  4 : mplp.openQ = atoi(optarg); break; */
        /* case 'p': mplp.flag |= MPLP_PER_SAMPLE; break; */

        default:
            fprintf(stderr, "Invalid option: '%c'\n", c);
            return 1;
        }
    }

    if (argc - optind != 3) return pileup_usage();
        
    char *samples_file = argv[optind];
    char *fasta_file = argv[optind + 1];
    char *pileup_file = argv[optind + 2];

    /* max_mem accounts for compressed bam plus pileup
       structures. max_input_mem accounts for the compressed bam
       size. */
    unsigned long max_input_mem = (float)opts.max_mem * 0.4;

    struct thread_queue *tq = 
        pileup_init(samples_file,
                    fasta_file,
                    pileup_file,
                    query_range_file,
                    opts.n_threads,
                    opts.n_readers,
                    max_input_mem,
                    opts.min_base_quality);

    printf("Starting input processing.\n");
    fflush(stdout);

    thread_queue_run(tq);
    thread_queue_free(tq);

    printf("Finished.\n");
    
    return 0;
}




/* provide thread_queue API functions */
void
pileup_worker(const struct managed_buf *in_bufs,
              unsigned more_input,
              struct managed_buf *out_bufs)
{
    struct managed_buf bam = { NULL, 0, 0 };
    struct managed_buf *out_buf = &out_bufs[0];
    unsigned s;
    for (s = 0; s != bam_samples.n; ++s) {
        bam_inflate(&in_bufs[s], &bam);
        pileup_tally_stats(bam, s);
    }

    if (! more_input)
        pileup_final_input();

    /* we need bqs and indel stats.  do not need basecall stats */
    for (s = 0; s != bam_samples.n; ++s) {
        pileup_prepare_bqs(s);
        pileup_prepare_indels(s);
    }

    if (bam.buf != NULL) free(bam.buf);

    struct pileup_data pdat = { 
        .calls = { NULL, 0, 0 },
        .quals = { NULL, 0, 0 },
        .n_match_lo_q = 0,
        .n_match_hi_q = 0,
        .n_indel = 0
    };

    struct pileup_locus_info ploc;
    char *out = out_buf->buf;
    out_buf->size = 0;
    while (pileup_next_pos()) {
        pileup_current_info(&ploc);
        for (s = 0; s != bam_samples.n; ++s) {
            pileup_current_data(s, &pdat);
            unsigned add = pdat.calls.size + pdat.quals.size + MAX_LABEL_LEN + 5;
            ALLOC_GROW_REMAP(out_buf->buf, out, out_buf->size + add, out_buf->alloc);
            out += 
                sprintf(out, "%s\t%s\t%u\t%c\t%u\t",
                        bam_samples.m[s].label,
                        ploc.refname,
                        ploc.pos,
                        ploc.refbase,
                        pdat.n_match_hi_q);
            
            strncpy(out, pdat.calls.buf, pdat.calls.size);
            out += pdat.calls.size;
            *out++ = '\t';
            strncpy(out, pdat.quals.buf, pdat.quals.size);
            out += pdat.quals.size;
            *out++ = '\n';
            out_buf->size = out - out_buf->buf;
        }
    }   

    /* frees statistics that have already been used in one of the
       distance calculations. */
    pileup_clear_finished_stats();
}


/* just write to stdout.  there is only one output file */
void
pileup_offload(void *par, const struct managed_buf *bufs)
{
    FILE *fh = par;
    fwrite(bufs[0].buf, 1, bufs[0].size, fh);
    fflush(fh);
}

void
pileup_on_create()
{
    batch_pileup_thread_init(bam_samples.n);
}


void
pileup_on_exit()
{
    batch_pileup_thread_free();
}


struct thread_queue *
pileup_init(const char *samples_file,
            const char *fasta_file,
            const char *pileup_file,
            const char *query_range_file,
            unsigned n_threads,
            unsigned n_readers,
            unsigned long max_input_mem,
            unsigned min_qual)
{
    bam_sample_info_init(samples_file, NULL);
    batch_pileup_init(min_qual, fasta_file);

    thread_params.pileup_fh = open_if_present(pileup_file, "w");

    unsigned long n_total_loci;
    thread_params.ranges = 
        parse_query_ranges(query_range_file,
                           &thread_params.n_ranges,
                           &n_total_loci);

    thread_params.reader_buf = malloc(n_readers * sizeof(struct bam_reader_par));
    thread_params.reader_par = malloc(n_readers * sizeof(void *));

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

    thread_queue_reader_t reader = { bam_reader, bam_scanner };
    struct thread_queue *tq =
        thread_queue_init(reader, thread_params.reader_par,
                          pileup_worker,
                          pileup_offload, thread_params.pileup_fh,
                          pileup_on_create,
                          pileup_on_exit,
                          n_threads, n_extra, n_readers, bam_samples.n,
                          n_outputs, max_input_mem);

    return tq;
}


void
pileup_free()
{
    bam_sample_info_free();
    batch_pileup_free();
    fclose(thread_params.pileup_fh);

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
