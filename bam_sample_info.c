/* Auxiliary helper functions for the dist program */

#include <pthread.h>

#include "bam_sample_info.h"
#include "common_tools.h"
#include "cache.h"
#include "khash.h"
#include "genome.h"

#include <assert.h>

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
    for (s = 0; s != bam_samples.n; ++s) {
        free(bam_samples.m[s].bam_file);
        fclose(bam_samples.m[s].fh);
    }
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
        ALLOC_GROW_TYPED(bam_sample_pairs.m, p + 1, alloc);
        bam_sample_pairs.m[p].s1 = kh_val(sample_map, itr[0]);
        bam_sample_pairs.m[p].s2 = kh_val(sample_map, itr[1]);
        ++p;
    }
    bam_sample_pairs.n = p;
    fclose(sample_pair_fh);
}
