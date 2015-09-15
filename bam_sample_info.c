/* Auxiliary helper functions for the dist program */

#include <pthread.h>

#include "bam_sample_info.h"
#include "common_tools.h"
#include "cache.h"
#include "khash.h"
#include "binomial_est.h"

#include <assert.h>

static pthread_mutex_t pair_stats_mtx = PTHREAD_MUTEX_INITIALIZER;

struct bam_sample_pair_info bam_sample_pairs;

struct bam_sample_info bam_samples;

KHASH_MAP_INIT_STR(remap_h, unsigned);


static void
init_sample_attributes(const char *samples_file, khash_t(remap_h) *sample_map);

static void
init_sample_pairs(const char *pair_file, khash_t(remap_h) *sample_map);

/* program-wide initialization of sample information */
void bam_sample_info_init(const char *samples_file, 
                          const char *sample_pairs_file)
{
    khash_t(remap_h) *sample_map = kh_init(remap_h);
    init_sample_attributes(samples_file, sample_map);
    init_sample_pairs(sample_pairs_file, sample_map);
    
    khiter_t k;
    for (k = kh_begin(sample_map); k != kh_end(sample_map); ++k)
        if (kh_exist(sample_map, k))
            free((char *)kh_key(sample_map, k));
    kh_destroy(remap_h, sample_map);
}

void bam_sample_info_free()
{
    free(bam_sample_pairs.m);
    unsigned s;
    for (s = 0; s != bam_samples.n; ++s)
        free(bam_samples.m[s].bam_file);
    free(bam_samples.m);

}

static void
init_sample_attributes(const char *samples_file, khash_t(remap_h) *sample_map)
{
    char label[1000], bam[1000];
    unsigned s = 0, n_fields, alloc = 0;
    khiter_t itr;
    int ret;
    
    FILE *samples_fh = open_if_present(samples_file, "r");
    while (! feof(samples_fh)) {
        n_fields = fscanf(samples_fh, "%s\t%s\n", label, bam);
        if (n_fields != 2) {
            fprintf(stderr, "Error, lines in samples file  %s must have two fields\n", 
                    samples_file);
            exit(1);
        }
        itr = kh_put(remap_h, sample_map, strdup(label), &ret);
        if (ret != 1 && ret != 2) {
            fprintf(stderr, "Found duplicate sample label %s while parsing samples file %s\n",
                    label, samples_file);
            exit(1);
        } else {
            kh_val(sample_map, itr) = s;
            ALLOC_GROW(bam_samples.m, s + 1, alloc);
            bam_samples.m[s].bam_file = strdup(bam);
            if (strlen(label) + 1 > sizeof(bam_samples.m[s].label)) {
                fprintf(stderr, "%s: error: sample label string must be "
                        "less than %Zu characters\n", __func__,
                        sizeof(bam_samples.m[s].label));
                exit(1);
            }
            strcpy(bam_samples.m[s].label, label);
            ++s;
        }
    }
    fclose(samples_fh);
    bam_samples.n = s;
}


/* parse a sample pairs file. */
static void
init_sample_pairs(const char *pair_file, khash_t(remap_h) *sample_map)
{
    int ret;
    FILE *sample_pair_fh = open_if_present(pair_file, "r");
    char label[2][MAX_LABEL_LEN];
    khiter_t itr[2];
    unsigned alloc = 0;

    if (! sample_pair_fh) {
        bam_sample_pairs.m = NULL;
        bam_sample_pairs.n = 0;
        return;
    }

    /* initialize the keyword REF as the reference sample */
    itr[0] = kh_put(remap_h, sample_map, strdup("REF"), &ret);
    assert(ret == 1 || ret == 2);
    kh_val(sample_map, itr[0]) = REFERENCE_SAMPLE;

    unsigned l, p = 0;
    while (! feof(sample_pair_fh)) {
        (void)fscanf(sample_pair_fh, "%s\t%s\n", label[0], label[1]);
        for (l = 0; l != 2; ++l)
            if ((itr[l] = kh_get(remap_h, sample_map, label[l])) == kh_end(sample_map)) {
                fprintf(stderr, "Error parsing %s file.  Couldn't find matching label %s in samples file\n",
                        pair_file, label[l]);
                exit(1);
            }
        ALLOC_GROW(bam_sample_pairs.m, p + 1, alloc);
        bam_sample_pairs.m[p].s1 = kh_val(sample_map, itr[0]);
        bam_sample_pairs.m[p].s2 = kh_val(sample_map, itr[1]);
        ++p;
    }
    bam_sample_pairs.n = p;
    fclose(sample_pair_fh);
}


/* increment all pair_dist_stats in bam_sample_pairs with stats */
void
accumulate_pair_stats(struct pair_dist_stats *stats)
{
    pthread_mutex_lock(&pair_stats_mtx);
    unsigned s, p;
    for (p = 0; p != bam_sample_pairs.n; ++p)
    {
        bam_sample_pairs.m[p].stats.total += stats[p].total;
        bam_sample_pairs.m[p].stats.confirmed_changed += stats[p].confirmed_changed;
        bam_sample_pairs.m[p].stats.cacheable += stats[p].cacheable;
        bam_sample_pairs.m[p].stats.cache_was_set += stats[p].cache_was_set;

        for (s = 0; s != N_STATS_CATEGORIES; ++s)
            bam_sample_pairs.m[p].stats.dist_count[s] += stats[p].dist_count[s];
    }

    pthread_mutex_unlock(&pair_stats_mtx);
}


void print_pair_stats(const char *stats_file)
{
    FILE *fh = open_if_present(stats_file, "w");
    pthread_mutex_lock(&pair_stats_mtx);
    fprintf(fh, 
            "%s\t%s\t%s\t%s\t%s\t%s", 
            "sample1", "sample2", "total", "cacheable",
            "cache_was_set", "confirmed_changed");

    unsigned s;
    for (s = 0; s != N_STATS_CATEGORIES; ++s)
        fprintf(fh, "\t%s", fuzzy_state_strings[s]);
    fprintf(fh, "\n");

    unsigned p;
    for (p = 0; p != bam_sample_pairs.n; ++p)
    {
        /* Print out statistics */
        unsigned
            s1 = bam_sample_pairs.m[p].s1,
            s2 = bam_sample_pairs.m[p].s2;
        fprintf(fh, "%s\t%s", 
                bam_samples.m[s1].label,
                s2 == REFERENCE_SAMPLE ? "REF" : bam_samples.m[s2].label);
        fprintf(fh, "\t%zu", bam_sample_pairs.m[p].stats.total);
        fprintf(fh, "\t%zu", bam_sample_pairs.m[p].stats.cacheable);
        fprintf(fh, "\t%zu", bam_sample_pairs.m[p].stats.cache_was_set);
        fprintf(fh, "\t%zu", bam_sample_pairs.m[p].stats.confirmed_changed);

        for (s = 0; s != N_STATS_CATEGORIES; ++s)
            fprintf(fh, "\t%zu", bam_sample_pairs.m[p].stats.dist_count[s]);
        fprintf(fh, "\n");
    }
    fclose(fh);
    pthread_mutex_unlock(&pair_stats_mtx);
}
