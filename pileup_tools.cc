#include "pileup_tools.h"
#include "tools.h"
#include "nucleotide_stats.h"
#include "file_utils.h"

#include <cassert>
#include <set>


char const PileupSummary::nucleotides[] = "ACGTN.acgtn,";

//maps every character to itself except the letters 'ACGTN.acgtn,'
//are mapped to 'N' (78)

char PileupSummary::code_to_redux[] = {
 '\x0', '\x1', '\x2', '\x3', '\x4', '\x5', '\x6', '\x7', '\x8', '\x9', '\xA', '\xB', '\xC', '\xD', '\xE', '\xF',
'\x10','\x11','\x12','\x13','\x14','\x15','\x16','\x17','\x18','\x19','\x1A','\x1B','\x1C','\x1D','\x1E','\x1F',
'\x20','\x21','\x22','\x23','\x24','\x25','\x26','\x27','\x28','\x29','\x2A','\x2B','\x4E','\x2D','\x4E','\x2F',
'\x30','\x31','\x32','\x33','\x34','\x35','\x36','\x37','\x38','\x39','\x3A','\x3B','\x3C','\x3D','\x3E','\x3F',
'\x40','\x4E','\x42','\x4E','\x44','\x45','\x46','\x4E','\x48','\x49','\x4A','\x4B','\x4C','\x4D','\x4E','\x4F',
'\x50','\x51','\x52','\x53','\x4E','\x55','\x56','\x57','\x58','\x59','\x5A','\x5B','\x5C','\x5D','\x5E','\x5F',
'\x60','\x4E','\x62','\x4E','\x64','\x65','\x66','\x4E','\x68','\x69','\x6A','\x6B','\x6C','\x6D','\x4E','\x6F',
'\x70','\x71','\x72','\x73','\x4E','\x75','\x76','\x77','\x78','\x79','\x7A','\x7B','\x7C','\x7D','\x7E','\x7F',
'\x80','\x81','\x82','\x83','\x84','\x85','\x86','\x87','\x88','\x89','\x8A','\x8B','\x8C','\x8D','\x8E','\x8F',
'\x90','\x91','\x92','\x93','\x94','\x95','\x96','\x97','\x98','\x99','\x9A','\x9B','\x9C','\x9D','\x9E','\x9F',
'\xA0','\xA1','\xA2','\xA3','\xA4','\xA5','\xA6','\xA7','\xA8','\xA9','\xAA','\xAB','\xAC','\xAD','\xAE','\xAF',
'\xB0','\xB1','\xB2','\xB3','\xB4','\xB5','\xB6','\xB7','\xB8','\xB9','\xBA','\xBB','\xBC','\xBD','\xBE','\xBF',
'\xC0','\xC1','\xC2','\xC3','\xC4','\xC5','\xC6','\xC7','\xC8','\xC9','\xCA','\xCB','\xCC','\xCD','\xCE','\xCF',
'\xD0','\xD1','\xD2','\xD3','\xD4','\xD5','\xD6','\xD7','\xD8','\xD9','\xDA','\xDB','\xDC','\xDD','\xDE','\xDF',
'\xE0','\xE1','\xE2','\xE3','\xE4','\xE5','\xE6','\xE7','\xE8','\xE9','\xEA','\xEB','\xEC','\xED','\xEE','\xEF',
'\xF0','\xF1','\xF2','\xF3','\xF4','\xF5','\xF6','\xF7','\xF8','\xF9','\xFA','\xFB','\xFC','\xFD','\xFE','\xFF'
};

int PileupSummary::quality_code_offset;

FastqType PileupSummary::fastq_type;

#define MIN(a,b) ((a) < (b) ? (a) : (b))

PileupSummary::PileupSummary(void) : 
    bases(NULL), bases_upper(NULL), bases_raw(NULL), quality_codes(NULL)
    {
        counts.raw_counts = NULL;
        counts.stats_index = NULL;
        counts.fbqs_cpd = NULL;
        counts.num_data = 0;
        memset(base_counts, 0, sizeof(base_counts[0]) * num_base_symbols);
        // memset(base_qual_sums, 0, sizeof(base_counts[0]) * num_base_symbols);
        sum_of_counts = 0;
    }


PileupSummary::~PileupSummary()
{
    if (bases) { delete[] bases; bases = NULL; }
    if (bases_upper) { delete[] bases_upper; bases_upper = NULL; }
    if (bases_raw) { delete[] bases_raw; bases_raw = NULL; }
    if (quality_codes) { delete[] quality_codes; quality_codes = NULL; }
    if (this->counts.raw_counts) { delete[] this->counts.raw_counts; this->counts.raw_counts = NULL; }
    if (this->counts.stats_index) { delete[] this->counts.stats_index; this->counts.stats_index = NULL; }
    if (this->counts.fbqs_cpd) { delete[] this->counts.fbqs_cpd; this->counts.fbqs_cpd = NULL; }

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
        sscanf(read_ptr, "%s\t%i\t%c\t%zi\t%n%*[^\t]\t%n", 
               this->reference, 
               &this->position, 
               &this->reference_base, 
               &span_depth, 
               &read_pos, &qual_pos);
    
    if (scanned_fields != 4)
    {
        fprintf(stderr, "PileupSummary::load_line: Warning: badly formatted line\n");
        assert(false);
        exit(10);
    }

    read_ptr += read_pos;

    if (this->quality_codes != NULL) delete[] this->quality_codes;
    this->quality_codes = new char[span_depth + 1];

    if (this->bases != NULL) delete[] this->bases;
    this->bases = new char[span_depth + 1];

    if (this->bases_upper != NULL) delete[] this->bases_upper;
    this->bases_upper = new char[span_depth + 1];

    if (this->bases_raw) delete[] this->bases_raw;
    
    int rawlen = qual_pos - read_pos - 1;
    this->bases_raw = new char[rawlen + 1]; // includes the extra space
    strncpy(this->bases_raw, read_ptr, rawlen);
    this->bases_raw[rawlen] = '\0';

    const int max_indel_size = 1000;
    char indel_sequence[max_indel_size + 1];

    char pileup_ccode;
    int indel_size;
    int indel_bin;
    int pileup_value;
    int current_read = 0;
    //int deletion_read = 0;

    while(*read_ptr != '\0' && *read_ptr != '\n')
    {
        pileup_ccode = *read_ptr;
        ++read_ptr;

        //reduce the pileup code
        int pileup_code = pileup_ccode;
        char pileup_redux = PileupSummary::code_to_redux[pileup_code];

        indel_size = 0;
        pileup_value = Nucleotide::base_to_index[pileup_code];
            
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
                /* 
                if (pileup_ccode == '.')
                {
                    real_base = toupper(this->reference_base);
                }
                else if (pileup_ccode == ',')
                {
                    real_base = tolower(this->reference_base);
                }
                else
                {
                    real_base = pileup_ccode;
                }
                */

                pileup_value = Nucleotide::base_to_index[static_cast<int>(real_base)];

                this->base_counts[pileup_value]++;
                this->sum_of_counts++;
                indel_bin = 0;

                this->bases_upper[current_read] = toupper(real_base);
                this->bases[current_read] = real_base;

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
        {
            //the end of a read.  ignored
        }
        else if (pileup_redux == '*')
        {
            /*
              see: http://samtools.sourceforge.net/samtools.shtml
              "Similarly, a pattern ‘-[0-9]+[ACGTNacgtn]+’ represents
              a deletion from the reference. The deleted bases will be
              presented as ‘*’ in the following lines."
             */
            this->bases[current_read++] = '*';
        }
        else if (pileup_redux == '\t' ||
                 pileup_redux == ' ')
        {
            //end of read_string, get the qual string
            sscanf(read_ptr, " %n", & read_pos); //eat white space
            read_ptr += read_pos;

            memcpy(this->quality_codes, read_ptr, span_depth);
            read_ptr += span_depth;

            this->quality_codes[span_depth] = '\0';

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
        if (this->bases[r] == '*')
        {
            ++num_seen_gaps;
            continue;
        }
        size_t write_pos = r - num_seen_gaps;
        this->bases[write_pos] = this->bases[r];
        this->bases_upper[write_pos] = this->bases_upper[r];
        this->quality_codes[write_pos] = this->quality_codes[r];
        // assert(this->quality_codes[r] != '~');
    }

    this->read_depth = span_depth;
    this->read_depth_match = span_depth - num_seen_gaps;
    this->bases[this->read_depth_match] = '\0';
    this->bases_upper[this->read_depth_match] = '\0';
    this->quality_codes[this->read_depth_match] = '\0';

}


size_t PileupSummary::quality(size_t read_num) const
{
    return QualityCodeToQuality(this->quality_codes[read_num],
                                PileupSummary::quality_code_offset);
}


// 1. encode the basecall, quality, strand combination as an integer
// initialize this->counts.raw_counts and this->counts.stats_index
void PileupSummary::parse(size_t min_quality_score)
{
    
    // populate raw_counts_flat with a sparse set of counts
    unsigned long *raw_counts_flat = new unsigned long[Nucleotide::num_bqs];
    unsigned long *rc_tmp = new unsigned long[Nucleotide::num_bqs];
    size_t *rc_ind_tmp = new size_t[Nucleotide::num_bqs];
    this->read_depth_high_qual = 0;

    std::fill(raw_counts_flat, raw_counts_flat + Nucleotide::num_bqs, 0);
    size_t rd = this->read_depth_match, nd = 0;

    for (size_t r = 0; r < rd; ++r)
    {
            
        size_t quality = this->quality(r);
        char basecall = this->bases[r];
        size_t basecall_index = Nucleotide::base_to_index[static_cast<size_t>(basecall)];
        if (basecall_index >= 4)
        {
            //since 'N' is common in pileup, but meaningless, we ignore it silently
            continue;
        }
        if (quality < min_quality_score)
        {
            continue;
        }
        size_t strand = isupper(basecall) ? Nucleotide::PLUS_STRAND : Nucleotide::MINUS_STRAND;
        size_t flat_index = Nucleotide::encode(basecall, quality, strand);
        raw_counts_flat[flat_index]++;
    }
    for (size_t flat_index = 0; flat_index != Nucleotide::num_bqs; ++flat_index)
    {
        if (raw_counts_flat[flat_index] != 0)
        {

            this->read_depth_high_qual += raw_counts_flat[flat_index];
            rc_tmp[nd] = raw_counts_flat[flat_index];
            rc_ind_tmp[nd] = flat_index;
            ++nd;
        }
    }

    this->counts.num_data = nd;
    if (this->counts.raw_counts != NULL)
        delete[] this->counts.raw_counts;

    this->counts.raw_counts = new unsigned long[nd];

    if (this->counts.stats_index != NULL)
        delete[] this->counts.stats_index;

    this->counts.stats_index = new size_t[nd];

    if (this->counts.fbqs_cpd != NULL)
        delete[] this->counts.fbqs_cpd;

    this->counts.fbqs_cpd = new double[nd * 4];

    std::copy(rc_tmp, rc_tmp + nd, this->counts.raw_counts);
    std::copy(rc_ind_tmp, rc_ind_tmp + nd, this->counts.stats_index);

    delete[] raw_counts_flat;
    delete[] rc_tmp;
    delete[] rc_ind_tmp;
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
            (*input->seen_codes_flags)[static_cast<size_t>(locus.quality_codes[rd])] = 1;
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
            max = locus.quality_codes[rd] > max ? locus.quality_codes[rd] : max;
            min = locus.quality_codes[rd] < min ? locus.quality_codes[rd] : min;
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


// determine the fastq offset of this pileup file
/*
FastqType FastqFileType(char const* pileup_file,
                        char * chunk_buffer_in,
                        size_t chunk_size,
                        size_t num_threads)
{
    FILE * pileup_input_fh = open_if_present(pileup_file, "r");
    CODES_FLAGS seen_codes_flags;
    std::fill(seen_codes_flags, seen_codes_flags + 256, 0);

    char * last_fragment;
    size_t max_pileup_line_size = 1000000;
    size_t bytes_wanted = chunk_size - max_pileup_line_size - 1;
    size_t bytes_read;
    fastq_type_input ** worker_input = new fastq_type_input *[num_threads];
    std::vector<char *> lines;
    CODES_FLAGS * seen_codes_flags_t = new CODES_FLAGS[num_threads];
    pthread_t * threads = new pthread_t[num_threads];

    for (size_t t = 0; t != num_threads; ++t)
    {
        std::fill(seen_codes_flags_t[t], seen_codes_flags_t[t] + 256, 0);
        worker_input[t] = 
            new fastq_type_input(t, lines.begin(), lines.end(), &seen_codes_flags_t[t]);
    }

    size_t fread_nsec;
    size_t total_bytes_read = 0;
    size_t total_fread_nsec = 0;

    while (! feof(pileup_input_fh))
    {
        bytes_read = FileUtils::read_until_newline(chunk_buffer_in, bytes_wanted,
                                                   max_pileup_line_size, pileup_input_fh,
                                                   &fread_nsec);

        lines = FileUtils::find_complete_lines_nullify(chunk_buffer_in, & last_fragment);
        total_bytes_read += bytes_read;
        total_fread_nsec += fread_nsec;

        for (size_t t = 0; t != num_threads; ++t)
        {
            size_t b = range_chunk_offset(0, lines.size(), num_threads, t, true);
            size_t e = range_chunk_offset(0, lines.size(), num_threads, t, false);

            worker_input[t]->beg = lines.begin() + b;
            worker_input[t]->end = lines.begin() + e;

            // fprintf(stdout, "thread %Zu: %Zu to %Zu of %Zu\n", t, b, e, lines.size());
            int rc = pthread_create(&threads[t], NULL, &fastq_type_worker, 
                                    static_cast<void *>(worker_input[t]));
            assert(rc == 0);
        }

        for (size_t t = 0; t != num_threads; ++t) {
            void *end;
            int rc = pthread_join(threads[t], &end);
            assert(0 == rc);
        }
        
    }
    fprintf(stderr, "FastqFileType: read %Zu bytes in %Zu ns.  (%5.3f MB/s)\n",
            total_bytes_read, total_fread_nsec,
            static_cast<float>(total_bytes_read) 
            / static_cast<float>(total_fread_nsec) / 1000.0);
            

    delete threads;

    fclose(pileup_input_fh);

    char seen_codes[256];
    char * seen_codes_ptr = seen_codes;
        
    for (size_t c = 0; c != 256; ++c)
    {
        for (size_t t = 0; t != num_threads; ++t)
        {
            if ((*worker_input[t]->seen_codes_flags)[c] == 1)
            {
                *seen_codes_ptr++ = static_cast<char>(c);
                break;
            }
        }
    }
    *seen_codes_ptr = '\0';

    for (size_t t = 0; t != num_threads; ++t)
    {
        delete worker_input[t];
    }
    delete worker_input;
    delete seen_codes_flags_t;

    fprintf(stderr, "Fastq codes seen in this file: %s\n", seen_codes);
    FastqType ftype = get_fastq_type(seen_codes);

    fprintf(stderr, "File type determined: %i\n", ftype);

    return ftype;
}
*/


void PileupSummary::SetFtype(FastqType _fastq_type) { 
    PileupSummary::fastq_type = _fastq_type; 
    PileupSummary::quality_code_offset = FastqTypeOffset(_fastq_type);
}


void PileupSummary::set_offset(int offset) {
    PileupSummary::quality_code_offset = offset;
}
