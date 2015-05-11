#include "pileup_tools.h"
#include "tools.h"
#include "nucleotide_stats.h"
#include "file_utils.h"
#include "defs.h"

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

static unsigned min_quality_score;

void pileup_init(unsigned _min_qual_score)
{
    if (_min_qual_score > 30)
    {
        fprintf(stderr, "Warning: minimum quality score of %u is probably too high...\n",
                _min_qual_score);
    }
    min_quality_score = _min_qual_score;
}

int PileupSummary::quality_code_offset;

FastqType PileupSummary::fastq_type;

#define MIN(a,b) ((a) < (b) ? (a) : (b))

PileupSummary::PileupSummary(void)
{
    bases.buf = bases_upper.buf = bases_raw.buf = quality_codes.buf = NULL;
    bases.alloc = bases_upper.alloc = bases_raw.alloc = quality_codes.alloc = 0;
    counts.num_data = 0;
    memset(base_counts_high_qual, 0, sizeof(base_counts_high_qual));
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
void PileupSummary::load_line(char *line)
{
    char *p = line;

    /* this->reference */
    unsigned ref_end = (p = (char *)memchr(p, '\t', MAX_PILEUP_LINE)) - line;
    strncpy(this->reference, line, ref_end);
    this->reference[ref_end] = '\0';

    /* this->position */
    this->position = strtol(++p, &p, 10);
    assert(*p++ == '\t');

    /* this->reference_base */
    this->reference_base = *p++;
    assert(*p++ == '\t');

    /* span_depth */
    unsigned span_depth = strtol(p, &p, 10);
    assert(*p++ == '\t');

    line = p;
    unsigned base_raw_len = (p = (char *)memchr(p, '\t', MAX_PILEUP_LINE)) - line;

    this->quality_codes.size = span_depth + 1;
    ALLOC_GROW_TYPED(this->quality_codes.buf, this->quality_codes.size,
                     this->quality_codes.alloc);

    this->bases.size = span_depth + 1;
    ALLOC_GROW_TYPED(this->bases.buf, this->bases.size,
                     this->bases.alloc);

    this->bases_upper.size = span_depth + 1;
    ALLOC_GROW_TYPED(this->bases_upper.buf, this->bases_upper.size,
                     this->bases_upper.alloc);

    this->bases_raw.size = base_raw_len + 1;
    ALLOC_GROW_TYPED(this->bases_raw.buf, this->bases_raw.size,
                     this->bases_raw.alloc);
    strncpy(this->bases_raw.buf, line, base_raw_len);
    this->bases_raw.buf[base_raw_len] = '\0';

    const int max_indel_size = 1000;
    char indel_sequence[max_indel_size + 1];

    char pileup_ccode;
    int indel_size, current_read = 0;
    //int deletion_read = 0;

    memset(this->base_counts_high_qual, 0, sizeof(this->base_counts_high_qual));

    this->counts.num_data = 0;
    this->insertions.clear();
    this->deletions.clear();

    while(*line != '\0' && *line != '\n')
    {
        pileup_ccode = *line++;

        //reduce the pileup code
        char pileup_redux = pileup_code_to_redux(pileup_ccode);

        indel_size = 0;

        // if (strchr("N+-", pileup_redux))
        if (pileup_redux == 'N' || pileup_redux == '+' || pileup_redux == '-')
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

                this->bases_upper.buf[current_read] = toupper(real_base);
                this->bases.buf[current_read] = real_base;

                ++current_read;
            }
            
            else 
            {
                //indel
                indel_size = strtol(line, &line, 10);

                // insertion
                memcpy(indel_sequence, line, MIN(indel_size, max_indel_size));
                indel_sequence[MIN(indel_size, max_indel_size)] = '\0';
                char *p = indel_sequence;
                while (*p != '\0') { *p = toupper(*p); ++p; }
                
                CHAR_MAP *container = pileup_redux == '+' ? &this->insertions : &this->deletions;
                CHAR_MAP::iterator it;
                ((it = container->find(indel_sequence)) == container->end())
                    ? (void)container->insert(CHAR_MAP::value_type(strdup(indel_sequence), 1))
                    : (void)(*it).second++;

                line += indel_size;
            }
        }

        else if (pileup_redux == '^')
        {
            //beginning of a read
            //ignore the next one character (a character-encoded mapping quality)
            // sscanf(line, "%*c");
            ++line;
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
            /* Not sure what I was thinking here...This seems unnecessary. */
            // sscanf(line, " %n", & read_pos); //eat white space
            // line += read_pos;

            memcpy(this->quality_codes.buf, line, span_depth);
            line += span_depth;

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
void PileupSummary::parse()
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
        this->base_counts_high_qual[basecall_index]++;
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
