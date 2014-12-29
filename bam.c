/* 
   magic char[4]
   n_ref int32_t
   n_bin int32_t
   bin uint32_t
   n_chunk int32_t
   chunk_beg uint64_t
   chunk_end uint64_t
   n_intv int32_t
   ioffset uint64_t
   n_no_coor uint64_t

*/

struct chunk_range_t {
    uint64_t chunk_beg, chunk_end;
};

struct bin_chunks_t {
    uint32_t bin;
    int32_t n_chunk;
    struct chunk_range_t *chunks;
};

struct contig_index_t {
    int32_t n_bin;
    struct bin_chunks_t *bins;
};


struct bam_index_t {
    char magic[4];
    int32_t n_ref;
    struct contig_index_t *contigs;
    int32_t n_intv;
    uint64_t *ioffset;
};
