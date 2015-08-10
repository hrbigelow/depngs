#ifndef _BAM_READER_H
#define _BAM_READER_H

#include "htslib/hts.h"
#include "htslib/bgzf.h"
#include "htslib/sam.h"

#include "cache.h"
#include "locus_range.h"

/* Provides read and scan functions for BAM files for use with
   thread_queue.  Enables user-specified locus ranges, and tandem
   reading of multiple BAM files. */

/* one instance of this per input sample per thread. */
struct bam_stats {
    hts_idx_t *idx;
    bam_hdr_t *hdr;
    BGZF *bgzf;
    hts_pair64_t *chunks;
    unsigned n_chunks;
    unsigned more_input; /* 1 means there is more input in this file. */
};

/* one instance per thread. parameter controlling scanning, and
   collecting certain information from the bam_scanner
   call. thread_queue passes this also to 'read' and 'worker'
   functions. */
struct bam_scanner_info {
    struct bam_stats *m;
    unsigned n;

    /* set of defined logical ranges */
    struct contig_region *qbeg, *qend;
    struct contig_span loaded_span;
};

/* called by up to n_readers threads at a time. par instructs the
   reader exactly what contents of each file to populate into
   bufs. returns 1 if there is more input available. */
unsigned
bam_reader(void *scanner_info, struct managed_buf *bufs);

/* called by at most one thread at a time.  reserves a logical range
   starting at internal position marker, and updates the marker to
   point to the end of the reserved range.  the range found will be
   the maximal range such that each sample's results are less than
   max_bytes.  updates par with settings to help accelerate the
   bam_reader call.  */
void
bam_scanner(void *scanner_info, size_t max_bytes);



/* return 1 if b overlaps any of the regions in [beg, end), 0
   otherwise. update (*beg) to point to first region that overlaps, or
   set (*beg) == end if no overlap found. */
int
rec_overlaps(bam1_t *b, 
             struct contig_region **beg,
             struct contig_region *end);


/* parse the next record of an uncompressed raw bam buffer into b,
   reallocating fields in b as necessary. return position of next bam
   record in the raw bam memory. adapted from
   htslib/sam.c:bam_read1. */
char *
bam_parse(char *bam_buf, bam1_t *b);


/* allocate a strndup'ed buffer of raw's contents of just one BAM
   record. (must be freed by the caller) */
char *
bam_duplicate_buf(char *raw);

/* inflate bgzf data stored in bgzf buffer, which contains the set of
   blocks defined by blocks and n_blocks.  store inflated BAM records
   in bam.  manage the size of bam.  only copy the portions of blocks
   defined by the virtual offsets. also, retrieves the serialized
   information about the actual loaded logical range into output
   argument loaded_range. */
void
bam_inflate(const struct managed_buf *bgzf,
            hts_pair64_t *chunks,
            unsigned n_chunks,
            struct managed_buf *bam);

/* initialize  */
void
bam_stats_init(const char *bam_file, struct bam_stats *bs);

/* release the resources of this bam_stat */
void
bam_stats_free(struct bam_stats *bs);


#endif /* _BAM_READER_H */
