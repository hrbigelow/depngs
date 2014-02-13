#include "pileup_tools.h"
#include "tools.h"
#include "nucleotide_stats.h"
#include "samutil/file_utils.h"

#include <cassert>
#include <set>


char const PileupSummary::nucleotides[] = "ACGTN.acgtn,";

//maps every character to itself except the letters 'ACGTN.acgtn,'
//are mapped to 'N' (78)
/*
char PileupSummary::code_to_redux[] = {
0, 1,  2,  3,  4,  5,  6,  7,  8,  9,  10,  11,  12,  13,  14,  15,
16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 
32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 78, 45, 78, 47, 
48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63, 
64, 78, 66, 78, 68, 69, 70, 78, 72, 73, 74, 75, 76, 77, 78, 79, 
80,81,82,83,  78,85,86,87,88,89,90,91,92,93,94,95, 
96, 78, 98, 78, 100,101,102, 78,104,105,106,107,108,109, 78,111, 
112,113,114,115, 78,117,118,119,120,121,122,123,124,125,126,127,
128,129,130,131,132,133,134,135,136,137,138,139,140,141,142,143,
144,145,146,147,148,149,150,151,152,153,154,155,156,157,158,159,
160,161,162,163,164,165,166,167,168,169,170,171,172,173,174,175,
176,177,178,179,180,181,182,183,184,185,186,187,188,189,190,191,
192,193,194,195,196,197,198,199,200,201,202,203,204,205,206,207,
208,209,210,211,212,213,214,215,216,217,218,219,220,221,222,223,
224,225,226,227,228,229,230,231,232,233,234,235,236,237,238,239,
240,241,242,243,244,245,246,247,248,249,250,251,252,253,254,255
};
*/

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


//load and parse a single line.  Assumes the line has all required fields.
//consumes the current line including newline character.
bool PileupSummary::load_line(char * read_pointer)
{

    size_t total_read_depth;
    int read_pos;

    int scanned_fields =
        sscanf(read_pointer, "%s\t%i\t%c\t%zi\t%n", this->_reference, 
               &this->_position, &this->_reference_base, 
               &total_read_depth, &read_pos);

    
    if (scanned_fields != 4)
    {
        fprintf(stderr, "PileupSummary::load_line: Warning: badly formatted line\n");
        return false;
    }

    read_pointer += read_pos;

    if (this->_quality_codes != NULL &&
        this->_bases != NULL &&
        this->_bases_upper != NULL)
    {
        delete this->_quality_codes;
        delete this->_bases;
        delete this->_bases_upper;
    }

    this->_quality_codes = new char[total_read_depth + 1];
    this->_bases = new char[total_read_depth + 1];
    this->_bases_upper = new char[total_read_depth + 1];

    const int max_indel_size = 1000;
    char indel_sequence[max_indel_size];

    char pileup_ccode;
    int indel_size;
    int indel_bin;
    int pileup_value;
    int current_read = 0;
    //int deletion_read = 0;

    while(*read_pointer != '\0')
    {
        pileup_ccode = *read_pointer;
        ++read_pointer;

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
                if (pileup_ccode == '.')
                {
                    real_base = toupper(this->_reference_base);
                }
                else if (pileup_ccode == ',')
                {
                    real_base = tolower(this->_reference_base);
                }
                else
                {
                    real_base = pileup_ccode;
                }

                pileup_value = Nucleotide::base_to_index[static_cast<int>(real_base)];

                this->base_counts[pileup_value]++;
                this->sum_of_counts++;
                indel_bin = 0;

                this->_bases_upper[current_read] = toupper(real_base);
                this->_bases[current_read] = real_base;

                ++current_read;
            }
            
            else 
            {
                //indel
                sscanf(read_pointer, "%i%n", &indel_size, &read_pos);
                read_pointer += read_pos;

                memcpy(indel_sequence, read_pointer, indel_size);
                read_pointer += indel_size;

                indel_sequence[indel_size] = '\0';

                // fgets(indel_sequence, indel_size + 1, file);
                for (int i = 0; i != indel_size; ++i)
                {
                    indel_sequence[i] = toupper(indel_sequence[i]);
                }
                indel_bin = std::min(this->_indel_histo_size, indel_size) * 
                    (pileup_redux == '+' ? 1 : -1);

                std::map<std::string, int> & indel = this->indel_seqs[indel_bin];
                std::string indel_str(indel_sequence);
                if (indel.find(indel_str) == indel.end())
                {
                    indel.insert(std::make_pair(indel_str, 0));
                }
                indel[indel_str]++;
            }

            //indel size histogram considers matches to be 'zero-length indels'
            this->indel_counts[indel_bin]++;
                
        }

        else if (pileup_redux == '^')
        {
            //beginning of a read
            //ignore the next one character (a character-encoded mapping quality)
            // sscanf(read_pointer, "%*c");
            ++read_pointer;
        }
        else if (pileup_redux == '$')
        {
            //the end of a read.  ignored
        }
        else if (pileup_redux == '*')
        {
            //what is this?  '*' is not documented, but according to Heng Li, represents a deletion.
            //fprintf(stderr, "asterisk found!\n");
            this->_bases[current_read++] = '*';
        }
        else if (pileup_redux == '\t' ||
                 pileup_redux == ' ')
        {
            //end of read_string, get the qual string
            sscanf(read_pointer, " %n", & read_pos); //eat white space
            read_pointer += read_pos;

            memcpy(this->_quality_codes, read_pointer, total_read_depth);
            read_pointer += total_read_depth;

            this->_quality_codes[total_read_depth] = '\0';
            // fgets(this->_quality_codes, total_read_depth + 1, file);

        }
        else 
        {
            fprintf(stderr, "Don't know this samtools pileup character code: %c\n", 
                    pileup_redux);
            return false;
        }
            
    }

    //now iterate through the bases and qualities, packing them
    size_t num_seen_gaps = 0;
    for (size_t read_pos = 0; read_pos != total_read_depth; ++read_pos)
    {
        if (this->_bases[read_pos] == '*')
        {
            ++num_seen_gaps;
            continue;
        }
        size_t write_pos = read_pos - num_seen_gaps;
        this->_bases[write_pos] = this->_bases[read_pos];
        this->_bases_upper[write_pos] = this->_bases_upper[read_pos];
        this->_quality_codes[write_pos] = this->_quality_codes[read_pos];
        assert(this->_quality_codes[read_pos] != '~');
    }

    this->_read_depth = total_read_depth - num_seen_gaps;
    this->_bases[this->_read_depth] = '\0';
    this->_bases_upper[this->_read_depth] = '\0';
    this->_quality_codes[this->_read_depth] = '\0';

    return true;
    
}


size_t PileupSummary::quality(size_t read_num) const
{
    return QualityCodeToQuality(this->_quality_codes[read_num],
                                PileupSummary::quality_code_offset);
}


FastqType PileupSummary::FastqFileType(char const* pileup_file,
                                       char * chunk_buffer_in,
                                       size_t chunk_size)
{
    FILE * pileup_input_fh = open_if_present(pileup_file, "r");
    char seen_codes_flags[256];
    std::fill(seen_codes_flags, seen_codes_flags + 256, 0);

    size_t nbytes_read, nbytes_unused = 0;
    char * last_fragment;
    char * read_pointer = chunk_buffer_in;

    while (! feof(pileup_input_fh))
    {
        nbytes_read = fread(read_pointer, 1, chunk_size - nbytes_unused, pileup_input_fh);
        
        std::vector<char *> lines =
            FileUtils::find_complete_lines_nullify(chunk_buffer_in, & last_fragment);
        
        read_pointer[nbytes_read] = '\0';
        
        for (size_t l = 0; l != lines.size(); ++l)
        {
            this->load_line(lines[l]);
            for (size_t rd = 0; rd != this->_read_depth; ++rd)
            {
                seen_codes_flags[static_cast<size_t>(this->_quality_codes[rd])] = 1;
            }
        }
        nbytes_unused = strlen(last_fragment);
        memmove(chunk_buffer_in, last_fragment, nbytes_unused);
        read_pointer = chunk_buffer_in + nbytes_unused;
    }
    fclose(pileup_input_fh);

    char seen_codes[256];
    char * seen_codes_ptr = seen_codes;
        
    for (size_t c = 0; c != 256; ++c)
    {
        if (seen_codes_flags[c] == 1)
        {
            *seen_codes_ptr++ = static_cast<char>(c);
        }
    }
    *seen_codes_ptr = '\0';

    fprintf(stderr, "Fastq codes seen in this file: %s\n", seen_codes);
    FastqType ftype = get_fastq_type(seen_codes);

    fprintf(stderr, "File type determined: %i\n", ftype);

    return ftype;
}


void PileupSummary::SetFtype(FastqType _fastq_type) { 
    PileupSummary::fastq_type = _fastq_type; 
    PileupSummary::quality_code_offset = FastqTypeOffset(_fastq_type);
}
