/* Convert multiple BAM files into a histogram of pileup statistics.
   Each locus is summarized by the numbers of majority and minority
   basecalls in each of several quality score bins.  This is called
   the summary statistic for the locus.  The number of loci for each
   distinct summary statistic is then reported. */
#include "bam_sample_info.h"
#include <stdlib.h>
#include <assert.h>
#include <pthread.h>
#include "htslib/sam.h"
#include "khash.h"
#include "thread_queue.h"
#include "chunk_strategy.h"
#include "batch_pileup.h"
#include "bam_reader.h"
#include "common_tools.h"
#include "locus_range.h"
#include "fasta.h"
#include "timer.h"

static struct {
    struct bam_filter_params bam_filter;
    struct batch_pileup_params bp_par;
    unsigned long max_mem;
    unsigned n_threads;
    unsigned n_max_reading;
    unsigned phred_offset;
} opts = { 
    { NULL, 5, 10, 0, 0 }, 
    {
        .skip_empty_loci = 0,
        .pseudo_depth = 1e6,
        .min_clash_qual = 20
    },
    .max_mem = 1e9,
    .n_threads = 1, 
    .n_max_reading = 1, 
    .phred_offset = 33
};


/* a place to hold allocated input for thread_queue that thread_queue
   does not own. */
static struct {
    struct bam_scanner_info *scanner_info_buf;
    void **reader_pars;
    unsigned n_threads;
    unsigned n_max_reading;
    struct contig_region *ranges;
    unsigned n_ranges;
    const char *fasta_file;
    unsigned pseudo_depth;
} thread_params;



int
pstats_usage()
{
    char *tmp_require = bam_flag2str(opts.bam_filter.rflag_require);
    char *tmp_filter  = bam_flag2str(opts.bam_filter.rflag_filter);

    fprintf(stderr,
            "\nUsage: dep pstats [options] samples.rdb in.fasta pstats.rdb\n\n"
            "Options:\n"
            "-t  INT    number of threads [%d]\n"
            "-R  INT    number of concurrent readers to allow (<= -t option) [%d]\n"
            "-m  INT    number bytes of memory to use [%ld]\n"
            "-F  INT    Output Phred Quality Encoding offset (33 or 64) [%d]\n"
            /* "-A  FLAG   do not discard anomalous read pairs\n" */
            /* "-B  FLAG   disable BAQ (per-Base Alignment Quality)\n" */
            /* "-C  INT    adjust mapping quality; recommended:50, disable:0 [0]\n" */
            /* "-E  FLAG   recalculate BAQ on the fly, ignore existing BQs\n" */
            /* "-G  STR    file list of read groups to exclude\n" */
            "-l  STR    file with locus ranges (chr start end) to process [empty = process all input]\n"
            "-q  INT    min mapping quality to include alignments [%d]\n"
            "-Q  INT    min base_quality/BAQ to include in pileup [%d]\n"
            "-x  STR    flags that must be present for inclusion (string or int) [%s]\n"
            "-X  STR    flags that must be absent for inclusion (string or int) [%s]\n"
            "\n\n"
            "samples.rdb has lines of <sample_label>\t</path/to/sample.bam>\n\n"
            "There must also exist an <in.fasta>.fai file produced by samtools faidx\n\n",
            opts.n_threads,
            opts.n_max_reading,
            opts.max_mem,
            opts.phred_offset,
            opts.bam_filter.min_map_quality,
            opts.bam_filter.min_base_quality,
            tmp_require,
            tmp_filter);
    
    free(tmp_require);
    free(tmp_filter);
    return 1;
}


/* max_input_mem is passed to thread_queue.  it is the full allotment
   of memory for all input buffers (n_samples * n_threads).  */
struct thread_queue *
pstats_init(const char *samples_file,
            const char *fasta_file,
            const char *query_range_file,
            unsigned n_threads,
            unsigned n_max_reading,
            unsigned long max_input_mem,
            struct bam_filter_params bam_filter,
            struct batch_pileup_params bp_par);

extern char *optarg;
extern int optind;


#define NBINS 9
#define QBIN_MAX 44

/* converts quality score to quality score bin */
static uint8_t q_to_qbin[] = { 0, 0, 0, 0, 0,
                               1, 1, 1, 1, 1,
                               2, 2, 2, 2, 2,
                               3, 3, 3, 3, 3,
                               4, 4, 4, 4, 4,
                               5, 5, 5, 5, 5,
                               6, 6, 6, 6, 6,
                               7, 7, 7, 7, 7,
                               8, 8, 8, 8, 8
};


typedef struct binned_qual_counts {
    struct {
        uint16_t major, minor;
    } bins[NBINS];
} binned_qual_counts_t;


#define MIN(a,b) ((a) < (b) ? (a) : (b))

khint32_t
binned_qual_hash_func(binned_qual_counts_t bqc)
{
    union {
        uint8_t maj[8];
        khint64_t key;
    } conv;
    unsigned i, lim = MIN(8, NBINS);
    for (i = 0; i != lim; ++i)
        conv.maj[i] = (uint8_t)bqc.bins[i].major;
    
    return kh_int64_hash_func(conv.key);
}


int
binned_qual_hash_equal(binned_qual_counts_t bqc1,
                       binned_qual_counts_t bqc2)
{
    return
        bqc1.bins[0].major == bqc2.bins[0].major
        && bqc1.bins[0].minor == bqc2.bins[0].minor
        && bqc1.bins[1].major == bqc2.bins[1].major
        && bqc1.bins[1].minor == bqc2.bins[1].minor
        && bqc1.bins[2].major == bqc2.bins[2].major
        && bqc1.bins[2].minor == bqc2.bins[2].minor
        && bqc1.bins[3].major == bqc2.bins[3].major
        && bqc1.bins[3].minor == bqc2.bins[3].minor
        && bqc1.bins[4].major == bqc2.bins[4].major
        && bqc1.bins[4].minor == bqc2.bins[4].minor
        && bqc1.bins[5].major == bqc2.bins[5].major
        && bqc1.bins[5].minor == bqc2.bins[5].minor
        && bqc1.bins[6].major == bqc2.bins[6].major
        && bqc1.bins[6].minor == bqc2.bins[6].minor
        && bqc1.bins[7].major == bqc2.bins[7].major
        && bqc1.bins[7].minor == bqc2.bins[7].minor
        && bqc1.bins[8].major == bqc2.bins[8].major
        && bqc1.bins[8].minor == bqc2.bins[8].minor;
    // return memcmp(&bqc1, &bqc2, sizeof(bqc1));
}


KHASH_INIT(bq_h, binned_qual_counts_t, khint64_t, 1, binned_qual_hash_func,
           binned_qual_hash_equal);

/* uncomment this at the end. messes up indentation. */
typedef khash_t(bq_h) bq_h_t;
                           
static __thread bq_h_t *tls_hash;
static pthread_mutex_t g_hash_mtx = PTHREAD_MUTEX_INITIALIZER;
static unsigned g_thread_num = 0;
static bq_h_t **g_part_hashes;

/* See hts.h.  This is the missing conversion table to go with
   seq_nt16_{table,str,int} */
const int seq_nt4_nt16[] = { 1, 2, 4, 8 };

/* convert *bqscount to the binned, base-ranked version. */
binned_qual_counts_t
bqs_count_to_key(struct bqs_count *bqs, unsigned n_bqs,
                 struct base_count bc)
{
    /* find the index of the majority base. */
    unsigned char mbi = 0;
    unsigned i, max_ct = 0;
    for (i = 0; i != 4; ++i)
        if (max_ct < bc.ct_filt[i]) {
            max_ct = bc.ct_filt[i];
            mbi = i;
        }
    /* The nt16 code for the major base */
    uint32_t nt16_maj = seq_nt4_nt16[mbi];
    unsigned qual_trunc, qbin;
    binned_qual_counts_t bqc;
    memset(bqc.bins, 0, sizeof(bqc.bins));
    
    for (i = 0; i != n_bqs; ++i) {
        qual_trunc = MIN(QBIN_MAX, bqs[i].qual);
        qbin = q_to_qbin[qual_trunc];
        if (bqs[i].base == nt16_maj)
            bqc.bins[qbin].major += bqs[i].ct;
        else
            bqc.bins[qbin].minor += bqs[i].ct;
    }
    return bqc;
}


/* must lock the mutex for global before calling this function. */
void
merge_histogram(bq_h_t *global, bq_h_t *part)
{
    khiter_t g_itr, p_itr;
    int empty;
    for (p_itr = kh_begin(part); p_itr != kh_end(part); ++p_itr) {
        if (! kh_exist(part, p_itr)) continue;
        g_itr = kh_put(bq_h, global, kh_key(part, p_itr), &empty);
        if (empty == 1)
            kh_val(global, g_itr) = kh_val(part, p_itr);
        else
            kh_val(global, g_itr) += kh_val(part, p_itr);
    }
}


struct histo_array {
    bq_h_t **parts;
    unsigned n_parts;
};


/* merge all histograms recursively into parts.  recursion rules:
   n_threads == 1: return
   n_threads == 2: 
 */
void *
recursive_merge_histogram(void *par)
{
    struct histo_array *mpar = par;
    unsigned n_threads = mpar->n_parts;
    bq_h_t **hashes = mpar->parts;

    if (n_threads <= 1)
        return NULL;
    if (n_threads == 2) {
        merge_histogram(hashes[0], hashes[1]);
        kh_destroy(bq_h, hashes[1]);
        return NULL;
    }
    
    pthread_t threads[2];
    struct histo_array parts[] = { 
        { hashes, n_threads / 2 },
        { hashes + n_threads / 2, n_threads - n_threads / 2 }
    };
    
    pthread_create(&threads[0], NULL, recursive_merge_histogram, &parts[0]);
    pthread_create(&threads[1], NULL, recursive_merge_histogram, &parts[1]);
    pthread_join(threads[0], NULL);
    pthread_join(threads[1], NULL);

    return NULL;
}


int main_pstats(int argc, char **argv)
{
    int c;
    char *query_range_file = NULL;
    int n_max_reading_set = 0;
    timer_init();
    
    /* adapted from samtools/bam_plcmd.c.*/
    while ((c = getopt(argc, argv, "t:R:m:F:ABC:EG:l:q:Q:x:X:")) >= 0) {
        switch (c) {
        case 't': opts.n_threads = strtol_errmsg(optarg, "-t (n_threads)"); break;
        case 'm': opts.max_mem = (size_t)strtod_errmsg(optarg, "-m (max_mem)"); break; 
        case 'R': 
            opts.n_max_reading = strtol_errmsg(optarg, "-R (n_max_reading)"); 
            n_max_reading_set = 1;
            break;
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
        case 'q': opts.bam_filter.min_map_quality = strtol_errmsg(optarg, "-q (min_map_quality)"); break;
        case 'Q': opts.bam_filter.min_base_quality = strtol_errmsg(optarg, "-Q (min_base_quality)"); break;

        case 'x' :
            opts.bam_filter.rflag_require = bam_str2flag(optarg);
            if ( opts.bam_filter.rflag_require<0 ) { fprintf(stderr,"Could not parse -x %s\n", optarg); return 1; }
            break;
        case 'X' :
            opts.bam_filter.rflag_filter = bam_str2flag(optarg);
            if ( opts.bam_filter.rflag_filter<0 ) { fprintf(stderr,"Could not parse --X %s\n", optarg); return 1; }
            break;

        /* case  4 : mplp.openQ = atoi(optarg); break; */
        /* case 'p': mplp.flag |= MPLP_PER_SAMPLE; break; */

        default:
            fprintf(stderr, "Invalid option: '%c'\n", c);
            return 1;
        }
    }

    /* by default, allow each thread to read unrestricted */
    if (! n_max_reading_set)
        opts.n_max_reading = opts.n_threads;

    if (argc - optind != 3) return pstats_usage();
        
    char *samples_file = argv[optind];
    char *fasta_file = argv[optind + 1];
    char *stats_file = argv[optind + 2];

    FILE *stats_fh = fopen(stats_file, "w");
    if (! stats_fh) {
        fprintf(stderr, "Error: Couldn't open statistics file %s for writing.\n",
                stats_file);
        exit(1);
    }

    /* max_mem accounts for compressed bam plus pileup structures plus
       output buffers.  max_input_mem accounts for the compressed bam
       size only. uncompressed it is about 3-4 times that size.
       printed out into pileup buffers it is about 4 times that size. */
    unsigned long max_input_mem = (float)opts.max_mem * 0.1;

    /* initialize our main structure */
    /* g_binqual_hash = kh_init(bq_h); */
    
    struct thread_queue *tq = 
        pstats_init(samples_file,
                    fasta_file,
                    stats_file,
                    query_range_file,
                    opts.n_threads,
                    opts.n_max_reading,
                    max_input_mem,
                    opts.bam_filter,
                    opts.bp_par);

    printf("Starting input processing.\n");
    fflush(stdout);

    thread_queue_run(tq);
    thread_queue_free(tq);

    /* This destroys all but g_part_hashes[0] */
    struct histo_array ha = { g_part_hashes, opts.n_threads };
    recursive_merge_histogram(&ha);

    unsigned i;
    khiter_t itr;
    
    /* hash has been merged. print it out. */
    for (itr = kh_begin(g_part_hashes[0]); itr != kh_end(g_part_hashes[0]); ++itr)
        if (kh_exist(g_part_hashes[0], itr)) {
            struct binned_qual_counts bqc = kh_key(g_part_hashes[0], itr);
            fprintf(stats_fh, "%Zi", kh_val(g_part_hashes[0], itr));
            for (i = 0; i != NBINS; ++i)
                fprintf(stats_fh, "\t%i", bqc.bins[i].major);
            for (i = 0; i != NBINS; ++i)
                fprintf(stats_fh, "\t%i", bqc.bins[i].minor);

            fprintf(stats_fh, "\n");
        }
    fclose(stats_fh);
    
    kh_destroy(bq_h, g_part_hashes[0]);
    
    printf("Finished.\n");
    
    return 0;
}

                                   
                                   
/* provide thread_queue API functions.
   vsi: cast this to struct bam_scanner_info. */
void
pstats_worker(struct managed_buf *in_bufs,
              unsigned more_input,
              void *vsi,
              struct managed_buf *out_bufs)
{
    struct managed_buf bam = { NULL, 0, 0 };
    struct bam_scanner_info *bsi = vsi;

    pileup_load_refseq_ranges(bsi);

    unsigned s;
    for (s = 0; s != bam_samples.n; ++s) {
        struct bam_stats *bs = &bsi->m[s];
        bam_inflate(&in_bufs[s], bs->chunks, bs->n_chunks, &bam);
        pileup_tally_stats(bam, bsi, s);
        in_bufs[s].size = 0;
        ALLOC_SHRINK(in_bufs[s].buf, in_bufs[s].size, in_bufs[s].alloc);
    }
    if (bam.buf != NULL) free(bam.buf);

    if (! more_input)
        pileup_final_input();

    /* we need bqs and indel stats.  do not need basecall stats */
    for (s = 0; s != bam_samples.n; ++s) {
        pileup_prepare_bqs(s);
        pileup_prepare_basecalls(s);
    }

    unsigned n_bqs_ct;
    struct bqs_count *bqs_ct = NULL;
    struct base_count bc;
    binned_qual_counts_t bqc;
    khiter_t itr;
    int empty;
    
    /* NOTE: this loop might have zero iterations if there is no data
       that intersected the region of interest.  */
    while (pileup_next_pos()) {
        for (s = 0; s != bam_samples.n; ++s) {
            pileup_current_bqs(s, &bqs_ct, &n_bqs_ct);
            bc = pileup_current_basecalls(s);
            bqc = bqs_count_to_key(bqs_ct, n_bqs_ct, bc);

            /* increment hash */
            itr = kh_put(bq_h, tls_hash, bqc, &empty);
            if (empty == 1)
                kh_val(tls_hash, itr) = 1;
            else
                ++kh_val(tls_hash, itr);
        }
    }   
    if (bqs_ct != NULL) free(bqs_ct);

    /* frees statistics that have already been used in one of the
       distance calculations. */
    pileup_clear_stats();
}


/* no operation done here since we are accumulating a global histogram
   as we go. */
void
pstats_offload(void *par, const struct managed_buf *bufs)
{
}

void
pstats_on_create()
{
    batch_pileup_thread_init(bam_samples.n, 
                             thread_params.fasta_file);
    pthread_mutex_lock(&g_hash_mtx);
    tls_hash = g_part_hashes[g_thread_num++];
    pthread_mutex_unlock(&g_hash_mtx);
}


void
pstats_on_exit()
{
    batch_pileup_thread_free();
    /* pthread_mutex_lock(&g_hash_mtx); */
    /* merge the hash */
    /* merge_histogram(g_binqual_hash, tls_hash); */
    /* pthread_mutex_unlock(&g_hash_mtx); */
    /* kh_destroy(bq_h, tls_hash); */
}

struct thread_queue *
pstats_init(const char *samples_file,
            const char *fasta_file,
            const char *locus_range_file,
            unsigned n_threads,
            unsigned n_max_reading,
            unsigned long max_input_mem,
            struct bam_filter_params bam_filter,
            struct batch_pileup_params bp_par)
{
    bam_sample_info_init(samples_file, NULL);
    bp_par.skip_empty_loci = 1;
    batch_pileup_init(bam_filter, bp_par);

    /* fasta_init is called from within batch_pileup_init */
    {
        unsigned t, s;
        const char **bam_files = malloc(bam_samples.n * sizeof(char *));
        for (s = 0; s != bam_samples.n; ++s)
            bam_files[s] = bam_samples.m[s].bam_file;
        
        struct bam_stats *all_stats =
            bam_stats_init_all(bam_files, bam_samples.n, n_threads);
        free(bam_files);
    
        thread_params.n_threads = n_threads;
        thread_params.scanner_info_buf = malloc(n_threads * sizeof(struct bam_scanner_info));
        thread_params.reader_pars = malloc(n_threads * sizeof(void *));

        for (t = 0; t != n_threads; ++t) {
            thread_params.scanner_info_buf[t] = (struct bam_scanner_info){
                .m = malloc(bam_samples.n * sizeof(struct bam_stats)),
                .n = bam_samples.n,
                .do_print_progress = 1
            };
            thread_params.reader_pars[t] = &thread_params.scanner_info_buf[t];
            for (s = 0; s != bam_samples.n; ++s)
                thread_params.scanner_info_buf[t].m[s] =
                    all_stats[t * bam_samples.n + s];
        }
        thread_params.fasta_file = fasta_file;
        free(all_stats);
    }
    
    /* estimated bytes left of bam input in a sample to switch to
       smaller chunking zones to avoid thread starvation.  see
       chunk_strategy.h */
    unsigned long bytes_zone2 = 1e8;
    unsigned long bytes_zone3 = 1e6;
    
    chunk_strategy_init(bam_samples.n, n_threads, 
                        locus_range_file, fasta_file,
                        bytes_zone2,
                        bytes_zone3);

    unsigned n_extra = 5;
    unsigned n_outputs = 0; /* output is indirectly produced thru
                               static hashes */

    thread_queue_reader_t reader = { bam_reader, bam_scanner };
    struct thread_queue *tq =
        thread_queue_init(reader,
                          thread_params.reader_pars,
                          pstats_worker,
                          pstats_offload, NULL,
                          pstats_on_create,
                          pstats_on_exit,
                          n_threads, n_extra, n_max_reading, bam_samples.n,
                          n_outputs, max_input_mem);

    g_part_hashes = malloc(sizeof(bq_h_t *) * n_threads);
    unsigned t;
    for (t = 0; t != n_threads; ++t)
        g_part_hashes[t] = kh_init(bq_h);

    return tq;
}


void
pstats_free()
{
    bam_sample_info_free();
    batch_pileup_free();

    unsigned t, s;
    for (t = 0; t != thread_params.n_threads; ++t) {
        for (s = 0; bam_samples.n; ++s)
            bam_stats_free(&thread_params.scanner_info_buf[t].m[s]);
        free(thread_params.scanner_info_buf[t].m);
    }

    free(thread_params.scanner_info_buf);
    free(thread_params.reader_pars);
    free(thread_params.ranges);
}
