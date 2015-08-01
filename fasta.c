#include "fasta.h"
#include "htslib/faidx.h"
#include "khash.h"

#include <stdio.h>
#include <assert.h>

static faidx_t *fasta_index; /* needed to efficiently retrieve sub-ranges of the
                         reference sequence. */

KHASH_MAP_INIT_STR(contig_h, unsigned);

/* does not own keys (keys owned by fasta_index) */
khash_t(contig_h) *contig_order;


/* initialize index (do not load any sequences).  call at start of
   program. */
void
fasta_init(const char *fasta_file)
{
    fasta_index = fai_load(fasta_file);
    if (fasta_index == NULL) {
        fprintf(stderr, "Error: %s:%u: Couldn't open fasta index file %s.fai\n",
                __FILE__, __LINE__, fasta_file);
        exit(1);
    }

    /* initialize contig_order hash */
    khiter_t itr;
    unsigned c, cn = faidx_nseq(fasta_index);
    int ret;
    const char *contig;
    contig_order = kh_init(contig_h);
    for (c = 0; c != cn; ++c) {
        contig = faidx_iseq(fasta_index, c);
        itr = kh_put(contig_h, contig_order, contig, &ret);
        kh_val(contig_order, itr) = c;
    }
}


/* free index.  call at end of program. */
void
fasta_free()
{
    kh_destroy(contig_h, contig_order);
}


/* return number of sequences in this fasta index */
int
fasta_nseq()
{
    return faidx_nseq(fasta_index);
}


/* return contig name from its tid */
const char *
fasta_get_contig(unsigned tid)
{
    return faidx_iseq(fasta_index, tid);
}


/* return tid from contig name, or -1 if not found. */
int
fasta_get_tid(const char *contig)
{
    khiter_t itr;
    if ((itr = kh_get(contig_h, contig_order, contig)) == kh_end(contig_order))
        return -1;
    else return kh_val(contig_order, itr);
}


int
fasta_seq_len(const char *contig)
{
    return faidx_seq_len(fasta_index, contig);
}


int
fasta_seq_ilen(unsigned tid)
{
    const char *contig = faidx_iseq(fasta_index, tid);
    return faidx_seq_len(fasta_index, contig);
}



/* fetch the sub-sequence of the fasta reference.  returns NULL on
   error. returned sequence must be freed by caller. */
char *
fasta_fetch_seq(const char *contig, int beg, int end)
{
    int fetch_len;
    char *seq =
        faidx_fetch_seq(fasta_index, contig, beg, end, &fetch_len);
    assert(fetch_len == end - beg);
    return seq;
}



/* fetch sub-sequence of the fasta reference, using tid to identify
   contig.  return NULL on error.  returned sequence must be freed by
   the caller. */
char *
fasta_fetch_iseq(unsigned tid, int beg, int end)
{
    int fetch_len;
    const char *contig = faidx_iseq(fasta_index, tid);
    char *seq = 
        faidx_fetch_seq(fasta_index, contig, beg, end, &fetch_len);

    if (! seq || fetch_len != (end - beg)) {
        fprintf(stderr, "Error: Couldn't retrieve reference subsequence %s:%u-%u\n",
                contig, beg, end);
        exit(1);
    }
    return seq;
}





