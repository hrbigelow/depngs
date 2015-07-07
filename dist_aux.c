/* Auxiliary helper functions for the dist program */

#include "dist_aux.h"
#include "common_tools.h"
#include "cache.h"
#include "khash.h"

#include <assert.h>

KHASH_MAP_INIT_STR(remap, int)

static khash_t(remap) *sample_map;
                      
void init_dist_aux()
{
    sample_map = kh_init(remap);
}

void free_dist_aux()
{
    khiter_t k;
    for (k = kh_begin(sample_map); k != kh_end(sample_map); ++k)
        if (kh_exist(sample_map, k)) free((char *)kh_key(sample_map, k));
    kh_destroy(remap, sample_map);
}

struct sample_attrs
init_sample_attributes(const char *samples_file)
{
    char jpd[1000], label[1000], pileup[1000];
    struct sample_attrs attrs;
    attrs.n = kh_size(sample_map);
    unsigned s, alloc = 0;
    ALLOC_GROW_TYPED(attrs.atts, attrs.n, alloc);
    khiter_t k;
    FILE *samples_fh = open_if_present(samples_file, "r");
    while (! feof(samples_fh))
    {
        (void)fscanf(samples_fh, "%s\t%s\t%s\n", label, jpd, pileup);
        if ((k = kh_get(remap, sample_map, label)) != kh_end(sample_map))
        {
            s = kh_val(sample_map, k);
            attrs.atts[s].file = strdup(pileup);
            nucleotide_stats_initialize(jpd, &attrs.atts[s].nuc_stats);
            if (strlen(label) + 1 > sizeof(attrs.atts[s].label)) 
            {
                fprintf(stderr, "%s: error: sample label string must be "
                        "less than %Zu characters\n", __func__,
                        sizeof(attrs.atts[s].label));
                exit(1);
            }
            strcpy(attrs.atts[s].label, label);
            attrs.atts[s].fh = fopen(pileup, "r");
            if (! attrs.atts[s].fh)
            {
                fprintf(stderr, "%s: error: couldn't open pileup input file %s\n",
                        __func__, pileup);
                exit(1);
            }
        }
    }
    fclose(samples_fh);

    return attrs;
}


/* parse a sample pairs file. */
struct sample_pairs
init_sample_pairs(const char *pair_file)
{
    struct sample_pairs spairs;

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
        ALLOC_GROW_TYPED(spairs.p, p + 1, alloc);
        spairs.p[p].s1 = kh_val(sample_map, k1);
        spairs.p[p].s2 = kh_val(sample_map, k2);
        ++p;
    }
    spairs.n = p;
    fclose(sample_pair_fh);

    k1 = kh_get(remap, sample_map, pseudo_key);
    kh_del(remap, sample_map, k1);

    return spairs;
}
