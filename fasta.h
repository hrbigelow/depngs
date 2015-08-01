#ifndef _FASTA_H
#define _FASTA_H

/* wrapper around faidx API to hide variables */

/* initialize index (do not load any sequences).  call at start of
   program. */
void
fasta_init(const char *fasta_file);


/* free index.  call at end of program. */
void
fasta_free();


/* return number of sequences in this fasta index */
int
fasta_nseq();


/* return contig name from its tid */
const char *
fasta_get_contig(unsigned tid);


/* return tid from contig name, or -1 if not found. */
int
fasta_get_tid(const char *contig);

/* get contig length, lookup by contig name */
int
fasta_seq_len(const char *contig);

/* get contig length, lookup by tid */
int
fasta_seq_ilen(unsigned tid);


/* fetch the sub-sequence of the fasta reference.  returns NULL on
   error. returned sequence must be freed by caller. */
char *
fasta_fetch_seq(const char *contig, int beg, int end);

/* fetch sub-sequence of the fasta reference, using tid to identify
   contig.  return NULL on error.  returned sequence must be freed by
   the caller. */
char *
fasta_fetch_iseq(unsigned tid, int beg, int end);


#endif /* _FASTA_H */
