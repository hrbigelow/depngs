#include "pileup_tools.h"
#include "tools.h"
#include "nucleotide_stats.h"
#include "file_utils.h"

#include <cassert>
#include <set>

/* maps every character to itself except the letters 'ACGTN.acgtn,' */
inline char pileup_code_to_redux(char code)
{
    switch(code)
    {
    case 'A': case 'C': case 'G': case 'T': case 'N': case '.':
    case 'a': case 'c': case 'g': case 't': case 'n': case ',':
        return 'N'; break;
    default: return code; break;
    }
}

int PileupSummary::quality_code_offset;

FastqType PileupSummary::fastq_type;

#define MIN(a,b) ((a) < (b) ? (a) : (b))

PileupSummary::PileupSummary(void)
{
    bases.buf = bases_upper.buf = bases_raw.buf = quality_codes.buf = NULL;
    bases.alloc = bases_upper.alloc = bases_raw.alloc = quality_codes.alloc = 0;
    counts.num_data = 0;
    memset(base_counts, 0, sizeof(base_counts[0]) * NUM_BASE_SYMBOLS);
    // memset(base_qual_sums, 0, sizeof(base_counts[0]) * num_base_symbols);
    sum_of_counts = 0;
    insertions = CHAR_MAP();
    deletions = CHAR_MAP();
}

PileupSummary::~PileupSummary()
{
    if (bases.alloc) { free(bases.buf); bases.alloc = 0; }
    if (bases_upper.alloc) { free(bases_upper.buf); bases_upper.alloc = 0; }
    if (bases_raw.alloc) { free(bases_raw.buf); bases_raw.alloc = 0; }
    if (quality_codes.alloc) { free(quality_codes.buf); quality_codes.alloc = 0; }

    CHAR_MAP::iterator it;
    for (it = this->insertions.begin(); it != this->insertions.end(); ++it)
        free((void *)(*it).first); // need free here since this is strdup'ed

    for (it = this->deletions.begin(); it != this->deletions.end(); ++it)
        free((void *)(*it).first);

}


/* load and parse a single line.  Assumes the line has all required
   fields.  Reads the input until the first newline or nul.
*/
void PileupSummary::load_line(const char *read_ptr)
{

    size_t span_depth;
    int read_pos, qual_pos;

    int scanned_fields =
        sscanf(read_ptr, 
               "%s\t%i\t%c\t%zi\t%n%*[^\t]\t%n", 
               this->reference, 
               &this->position, 
               &this->reference_base, 
               &span_depth, 
               &read_pos, &qual_pos);
    
    if (scanned_fields != 4)
    {
        char *line = strndup(read_ptr, strchr(read_ptr, '\n') - read_ptr);
        fprintf(stderr, "PileupSummary::load_line: Warning: badly formatted line:\n%s\n",
                line);
        free(line);
        assert(false);
        exit(10);
    }

    read_ptr += read_pos;

    this->quality_codes.size = span_depth + 1;
    ALLOC_GROW_TYPED(this->quality_codes.buf, this->quality_codes.size,
                     this->quality_codes.alloc);

    this->bases.size = span_depth + 1;
    ALLOC_GROW_TYPED(this->bases.buf, this->bases.size,
                     this->bases.alloc);

    this->bases_upper.size = span_depth + 1;
    ALLOC_GROW_TYPED(this->bases_upper.buf, this->bases_upper.size,
                     this->bases_upper.alloc);

    int rawlen = qual_pos - read_pos - 1;
    this->bases_raw.size = rawlen + 1;
    ALLOC_GROW_TYPED(this->bases_raw.buf, this->bases_raw.size,
                     this->bases_raw.alloc);
    strncpy(this->bases_raw.buf, read_ptr, rawlen);
    this->bases_raw.buf[rawlen] = '\0';

    const int max_indel_size = 1000;
    char indel_sequence[max_indel_size + 1];

    char pileup_ccode;
    int indel_size, pileup_value, current_read = 0;
    //int deletion_read = 0;

    memset(this->base_counts, 0, sizeof(this->base_counts));
    this->sum_of_counts = 0;

    this->counts.num_data = 0;
    this->insertions.clear();
    this->deletions.clear();

    while(*read_ptr != '\0' && *read_ptr != '\n')
    {
        pileup_ccode = *read_ptr;
        ++read_ptr;

        //reduce the pileup code
        int pileup_code = pileup_ccode;
        char pileup_redux = pileup_code_to_redux(pileup_ccode);

        indel_size = 0;
        pileup_value = base_to_index(pileup_code);
            
        if (strchr("N+-", pileup_redux))
        // if (pileup_redux == 'N' || pileup_redux == '+' || pileup_redux == '-')
        {
            //actual sequence

            if (pileup_redux == 'N')
            { 
                // ACGTN.acgtn, are mapped to 'N'
                char real_base;
                switch(pileup_ccode) {
                case '.': real_base = toupper(this->reference_base); break;
                case ',': real_base = tolower(this->reference_base); break;
                default : real_base = pileup_ccode; break;
                }

                pileup_value = base_to_index(real_base);

                this->base_counts[pileup_value]++;
                this->sum_of_counts++;
                this->bases_upper.buf[current_read] = toupper(real_base);
                this->bases.buf[current_read] = real_base;

                ++current_read;
            }
            
            else 
            {
                //indel
                sscanf(read_ptr, "%i%n", &indel_size, &read_pos);
                read_ptr += read_pos;

                // insertion
                memcpy(indel_sequence, read_ptr, MIN(indel_size, max_indel_size));
                indel_sequence[MIN(indel_size, max_indel_size)] = '\0';
                char *p = indel_sequence;
                while (*p != '\0') { *p = toupper(*p); ++p; }
                
                CHAR_MAP *container = pileup_redux == '+' ? &this->insertions : &this->deletions;
                CHAR_MAP::iterator it;
                ((it = container->find(indel_sequence)) == container->end())
                    ? (void)container->insert(CHAR_MAP::value_type(strdup(indel_sequence), 1))
                    : (void)(*it).second++;

                read_ptr += indel_size;
            }
        }

        else if (pileup_redux == '^')
        {
            //beginning of a read
            //ignore the next one character (a character-encoded mapping quality)
            // sscanf(read_ptr, "%*c");
            ++read_ptr;
        }
        else if (pileup_redux == '$')
            ; /* the end of a read.  ignored */


        /* see: http://samtools.sourceforge.net/samtools.shtml
           "Similarly, a pattern ‘-[0-9]+[ACGTNacgtn]+’ represents a
           deletion from the reference. The deleted bases will be
           presented as ‘*’ in the following lines."
        */
        else if (pileup_redux == '*')
            this->bases.buf[current_read++] = '*';

        else if (pileup_redux == '\t' || pileup_redux == ' ')
        {
            //end of read_string, get the qual string
            sscanf(read_ptr, " %n", & read_pos); //eat white space
            read_ptr += read_pos;

            memcpy(this->quality_codes.buf, read_ptr, span_depth);
            read_ptr += span_depth;

            this->quality_codes.buf[span_depth] = '\0';

        }
        else 
        {
            fprintf(stderr, "Don't know this samtools pileup character code: %c\n", 
                    pileup_redux);
            exit(10);
        }
            
    }

    // now iterate through the bases and qualities, packing them to
    // only include those bases that are in a CIGAR 'M' state.
    size_t num_seen_gaps = 0;
    for (size_t r = 0; r != span_depth; ++r)
    {
        if (this->bases.buf[r] == '*')
        {
            ++num_seen_gaps;
            continue;
        }
        size_t write_pos = r - num_seen_gaps;
        this->bases.buf[write_pos] = this->bases.buf[r];
        this->bases_upper.buf[write_pos] = this->bases_upper.buf[r];
        this->quality_codes.buf[write_pos] = this->quality_codes.buf[r];
        // assert(this->quality_codes.buf[r] != '~');
    }

    this->read_depth = span_depth;
    this->read_depth_match = span_depth - num_seen_gaps;
    this->bases.buf[this->read_depth_match] = '\0';
    this->bases_upper.buf[this->read_depth_match] = '\0';
    this->quality_codes.buf[this->read_depth_match] = '\0';

}


size_t PileupSummary::quality(size_t read_num) const
{
    return QualityCodeToQuality(this->quality_codes.buf[read_num],
                                PileupSummary::quality_code_offset);
}


// 1. encode the basecall, quality, strand combination as an integer
// initialize this->counts.raw_counts and this->counts.stats_index
void PileupSummary::parse(size_t min_quality_score)
{
    // populate raw_counts_flat with a sparse set of counts
    unsigned long raw_counts_flat[NUC_NUM_BQS];
    this->read_depth_high_qual = 0;

    std::fill(raw_counts_flat, raw_counts_flat + NUC_NUM_BQS, 0);
    size_t rd = this->read_depth_match, nd = 0;
    size_t fi;

    for (size_t r = 0; r < rd; ++r)
    {
            
        size_t quality = this->quality(r);
        char basecall = this->bases.buf[r];
        size_t basecall_index = base_to_index(basecall);

        /* since 'N' is common in pileup, but meaningless, ignore
           silently */
        if (basecall_index >= 4) continue;
        if (quality < min_quality_score) continue;

        size_t strand = isupper(basecall) ? NUC_PLUS_STRAND : NUC_MINUS_STRAND;
        fi = encode_nucleotide(basecall, quality, strand);
        raw_counts_flat[fi]++;
    }
    for (fi = 0; fi != NUC_NUM_BQS; ++fi)
    {
        if (raw_counts_flat[fi] != 0)
        {
            this->read_depth_high_qual += raw_counts_flat[fi];
            this->counts.stats[nd].ct = raw_counts_flat[fi];
            this->counts.stats_index[nd] = fi;
            ++nd;
        }
    }
    this->counts.num_data = nd;
}


typedef char CODES_FLAGS[256];

struct fastq_type_input
{
    size_t thread_num;
    std::vector<char *>::iterator beg;
    std::vector<char *>::iterator end;
    CODES_FLAGS * seen_codes_flags; // do not own this
    fastq_type_input(size_t thread_num,
                     std::vector<char *>::iterator beg,
                     std::vector<char *>::iterator end,
                     CODES_FLAGS * seen_codes_flags) :
        thread_num(thread_num), beg(beg), end(end), seen_codes_flags(seen_codes_flags)
    {
    }
    
};

void * fastq_type_worker(void * args)
{
    fastq_type_input * input = static_cast<fastq_type_input *>(args);

    PileupSummary locus;
    std::vector<char *>::iterator it;
    for (it = input->beg; it != input->end; ++it)
    {
        locus.load_line(*it);
        for (size_t rd = 0; rd != locus.read_depth_match; ++rd)
        {
            (*input->seen_codes_flags)[static_cast<size_t>(locus.quality_codes.buf[rd])] = 1;
        }
    }
    pthread_exit((void*) 0);
}


// return the probable fastq offset (33 or 64) for this pileup file.
// scan up to chunk_size bytes.
// 
int fastq_offset(const char *pileup_file,
                 char *chunk_buffer_in,
                 size_t chunk_size)
{
    unsigned char min = 255, max = 0;
    FILE *pileup_input_fh = open_if_present(pileup_file, "r");
    if (! pileup_input_fh)
        return -1;

    size_t fread_nsec;

    FileUtils::read_until_newline(chunk_buffer_in, chunk_size, 1e6, 
                                  pileup_input_fh, &fread_nsec);
    
    char *last_fragment;
    std::vector<char *> lines = 
        FileUtils::find_complete_lines_nullify(chunk_buffer_in, &last_fragment);

    std::vector<char *>::iterator it;
    PileupSummary locus;
    for (it = lines.begin(); it != lines.end(); ++it)
    {
        locus.load_line(*it);
        for (size_t rd = 0; rd != locus.read_depth_match; ++rd)
        {
            max = locus.quality_codes.buf[rd] > max ? locus.quality_codes.buf[rd] : max;
            min = locus.quality_codes.buf[rd] < min ? locus.quality_codes.buf[rd] : min;
        }
    }
    fclose(pileup_input_fh);

    // reject probable Solexa format
    if (min >= 59 && min < 64)
    {
        fprintf(stderr, "Warning: quality code minimum is %i, which is likely the unsupported Solexa format\n",
                min);
        return -1;
    }
    // Sanger and Illumina18 format range [33, 126]
    else if (min < 64 && max <= 126)
        return 33;

    // Illumina13 and Illumina15  have range [64, 126]
    else if (min >= 64 && max <= 126)
        return 64;

    else
    {
        fprintf(stderr, "Warning: quality code range is from %i to %i, which doesn't fit any known encoding\n",
                min, max);
        return -1;
    }

}


void PileupSummary::SetFtype(FastqType _fastq_type) { 
    PileupSummary::fastq_type = _fastq_type; 
    PileupSummary::quality_code_offset = FastqTypeOffset(_fastq_type);
}


void PileupSummary::set_offset(int offset) {
    PileupSummary::quality_code_offset = offset;
}
