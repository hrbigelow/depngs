#ifndef _BAM_READER_H
#define _BAM_READER_H

#include "htslib/hts.h"
#include "htslib/bgzf.h"
#include "htslib/sam.h"

#include "cache.h"

/* Provides read and scan functions for BAM files for use with
   thread_queue.  Enables user-specified locus ranges, and tandem
   reading of multiple BAM files. */

/* One instance for each of thread_queue's n_readers
   thread_queue_reader_t instances. It retains state between calls of
   bam_reader and bam_scanner.  */
struct bam_reader_par {
    hts_idx_t *idx;
    unsigned n_idx;
    
    BGZF **bgzf;

    /* set of defined logical ranges */
    struct pair_ordering_range *qbeg, *qend;

    hts_pair64_t **chunks; /* chunks[s]: set of chunks for sample s */
    unsigned *n_chunks; /* n_chunks[s]: number of chunks for sample s */
};

/* called by up to n_readers threads at a time. par instructs the
   reader exactly what contents of each file to populate into bufs. */
void bam_reader(void *par, struct managed_buf *bufs);

/* called by at most one thread at a time.  reserves a logical range
   starting at internal position marker, and updates the marker to
   point to the end of the reserved range.  the range found will be
   the maximal range such that each sample's results are less than
   max_bytes.  updates par with settings to help accelerate the
   bam_reader call.  */
void bam_scanner(void *par, unsigned max_bytes);


/* parse the next record of an uncompressed raw bam buffer into b,
   reallocating fields in b as necessary. return position of next bam
   record in the raw bam memory. adapted from
   htslib/sam.c:bam_read1. */
char *bam_parse(char *bam_buf, bam1_t *b);

#endif /* _BAM_READER_H */
