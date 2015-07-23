/* Auxiliary helper functions for the dist program */

#include <pthread.h>

#include "bam_sample_info.h"
#include "common_tools.h"
#include "cache.h"
#include "khash.h"

#include <assert.h>

struct bam_sample_pair_info bam_sample_pairs;

struct bam_sample_info bam_samples;

KHASH_MAP_INIT_STR(remap, int);


static void
init_sample_attributes(const char *samples_file, khash_t(remap) *sample_map);

static void
init_sample_pairs(const char *pair_file, khash_t(remap) *sample_map);

/* program-wide initialization of sample information */
void bam_sample_info_init(const char *samples_file, 
                          const char *sample_pairs_file)
{
    khash_t(remap) *sample_map = kh_init(remap);
    init_sample_pairs(sample_pairs_file, sample_map);
    init_sample_attributes(samples_file, sample_map);
    
    khiter_t k;
    for (k = kh_begin(sample_map); k != kh_end(sample_map); ++k)
        if (kh_exist(sample_map, k)) free((char *)kh_key(sample_map, k));
    kh_destroy(remap, sample_map);
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
init_sample_attributes(const char *samples_file, khash_t(remap) *sample_map)
{
    char label[1000], bam[1000];
    bam_samples.n = kh_size(sample_map);
    unsigned s, alloc = 0;
    ALLOC_GROW_TYPED(bam_samples.m, bam_samples.n, alloc);
    khiter_t k;
    FILE *samples_fh = open_if_present(samples_file, "r");
    while (! feof(samples_fh)) {
        (void)fscanf(samples_fh, "%s\t%s\n", label, bam);
        if ((k = kh_get(remap, sample_map, label)) != kh_end(sample_map)) {
            s = kh_val(sample_map, k);
            bam_samples.m[s].bam_file = strdup(bam);
            if (strlen(label) + 1 > sizeof(bam_samples.m[s].label)) {
                fprintf(stderr, "%s: error: sample label string must be "
                        "less than %Zu characters\n", __func__,
                        sizeof(bam_samples.m[s].label));
                exit(1);
            }
            strcpy(bam_samples.m[s].label, label);
        }
    }
    fclose(samples_fh);
}


/* parse a sample pairs file. */
static void
init_sample_pairs(const char *pair_file, khash_t(remap) *sample_map)
{
    unsigned idx = 0;
    int ret;
    FILE *sample_pair_fh = open_if_present(pair_file, "r");
    char label[2][MAX_LABEL_LEN];
    khiter_t k1, k2;
    unsigned alloc = 0;

    /* initialize the keyword PSEUDO as the reference sample */
    char pseudo_key[] = "REF";
    k1 = kh_put(remap, sample_map, pseudo_key, &ret);
    assert(ret == 0 || ret == 1);
    kh_val(sample_map, k1) = -1;

    unsigned p = 0;
    while (! feof(sample_pair_fh))
    {
        (void)fscanf(sample_pair_fh, "%s\t%s\n", label[0], label[1]);
        if ((k1 = kh_get(remap, sample_map, label[0])) == kh_end(sample_map))
        {
            k1 = kh_put(remap, sample_map, strdup(label[0]), &ret);
            assert(ret == 0 || ret == 1);
            kh_val(sample_map, k1) = idx++;
        }
        if ((k2 = kh_get(remap, sample_map, label[1])) == kh_end(sample_map))
        {
            k2 = kh_put(remap, sample_map, strdup(label[1]), &ret);
            assert(ret == 0 || ret == 1);
            kh_val(sample_map, k2) = idx++;
        }
        ALLOC_GROW_TYPED(bam_sample_pairs.m, p + 1, alloc);
        bam_sample_pairs.m[p].s1 = kh_val(sample_map, k1);
        bam_sample_pairs.m[p].s2 = kh_val(sample_map, k2);
        ++p;
    }
    bam_sample_pairs.n = p;
    fclose(sample_pair_fh);

    k1 = kh_get(remap, sample_map, pseudo_key);
    kh_del(remap, sample_map, k1);
}
