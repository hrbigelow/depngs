#include "tools.h"

#include <cstring>
#include <algorithm>
#include <cmath>
#include <cassert>
#include <cstdarg>
#include <map>


//Create a lookup table to go from char to index, given an
//index-to-char map expressed in a char *
void CharLookupTable(char const* index_to_char_map, int default_value, int *lookup_table)
{
    //!!! this doesn't do what's expected
    std::fill(lookup_table, lookup_table + 256, default_value);
    for (int ci = 0; index_to_char_map[ci] != '\0'; ++ci)
    {
        lookup_table[static_cast<int>(index_to_char_map[ci])] = ci;
    }
}


    
//compute log(n1 + n2), where n1 or n2 are underflow numbers expressed
//in log base 2.  If the difference between log_n1 and log_n2 is
//greater than 32, the lesser number is treated as log2(0.0)
float SumLog2NumbersTruncated(float log_n1, float log_n2)
{
    
    if (fabsf(log_n1 - log_n2) > 32)
    {
        return std::max(log_n1, log_n2);
    }
    if (log_n1 == -INFINITY && log_n2 == -INFINITY)
    {
        return -INFINITY;
    }
    else 
    {
        float log_min = floorf(std::min(log_n1, log_n2));
        float sum_scaled = 
            (exp2f(log_n1 - log_min) + exp2f(log_n2 - log_min));

        return log2f(sum_scaled) + log_min;
    }
}


//the accumulator function for Log2Accumulate
std::pair<float, SAMPLE_POINT>
Log2AccumulateStep(float const& sum_log_numbers, SAMPLE_POINT const& input)
{
    float sum = SumLog2NumbersTruncated(sum_log_numbers, input.second);
    return std::make_pair(sum, std::make_pair(input.first, sum));
}


//produce a new set of sample points at the same X values, representing
//the cumulant of the sample points in [begin, end)
SAMPLE Log2Accumulate(SAMPLE::iterator begin, SAMPLE::iterator end)
{
    return MapFold<SAMPLE, SAMPLE, float>(-INFINITY, begin, end, 
                                          &Log2AccumulateStep);
}


//accumulator function for calculating an upper bound on the norm
//returns a pair of accumulator, output point
std::pair<float, SAMPLE_POINT>
UpperBoundStep(float const& sum_log_numbers, 
               SAMPLE_POINT const& point1,
               SAMPLE_POINT const& point2)
{
    float log_area = 
        std::max(point1.second, point2.second) +
        log2f(point2.first - point1.first);
        
    float sum = SumLog2NumbersTruncated(sum_log_numbers, log_area);
    return std::make_pair(sum, std::make_pair(point2.first, sum));
}


//produce a new set of sample points at the same X values, representing
//an upper bound on the Norm
SAMPLE UpperBoundNormLog2(SAMPLE::iterator begin, SAMPLE::iterator end)
{
    return MapFold2<SAMPLE, SAMPLE, float>(-INFINITY, begin, end, 
                                          &UpperBoundStep);
}



std::pair<float, SAMPLE_POINT>
LowerBoundLog2Step(float const& sum_log_numbers, 
                   SAMPLE_POINT const& point1,
                   SAMPLE_POINT const& point2)
{
    float log_area = 
        std::min(point1.second, point2.second) +
        log2f(point2.first - point1.first);
    
    float sum = SumLog2NumbersTruncated(sum_log_numbers, log_area);
    return std::make_pair(sum, std::make_pair(point2.first, sum));
}


//produce a new set of sample points at the same X values, representing
//the cumulative area of the lower bound step function.
//the provided curve, given by points [begin, end) is in log form
//the desired cumulative lower bound curve is also in log form,
//but represents the integral of the non-log form lower-bound step functoin
SAMPLE LowerBoundNormLog2(SAMPLE::iterator begin, SAMPLE::iterator end)
{
    return MapFold2<SAMPLE, SAMPLE, float>(-INFINITY, begin, end, 
                                           &LowerBoundLog2Step);
}


char QualityToQualityCode(int quality, size_t quality_code_offset)
{
    return static_cast<char>(quality + quality_code_offset);
}



size_t QualityCodeToQuality(char quality_code, size_t quality_code_offset)
{
    return static_cast<size_t>(quality_code) - quality_code_offset
;
}


size_t FastqTypeOffset(FastqType ftype)
{
    size_t offset = 0;
    switch(ftype)
    {
    case Sanger: case Illumina18: offset = 33; break;
    case Solexa: case Illumina13: case Illumina15: offset = 64; break;
    case None: offset = 0; break;
    }
    return offset;
}


//Translates Phred quality into error probability
float QualityToErrorProb(int quality)
{
    return exp10(- static_cast<float>(quality) / 10.0);
}


int error_prob_to_quality(float error_prob)
{
    // !!! safe?
    double quality_est = (log2f(error_prob) / log2f(10.0)) * -10.0;
    int quality_ceil = static_cast<int>(ceil(quality_est));
    int quality_floor = static_cast<int>(floor(quality_est));

    double ceil_diff = quality_ceil - quality_est;
    double floor_diff = quality_est - quality_floor;
    return ceil_diff < floor_diff ? quality_ceil : quality_floor;

}


//prints out a line in pileup format
//
void PrintPileupLine(FILE * out_fh,
                     char const* reference,
                     int position,
                     char reference_base,
                     char const* start_base,
                     int const* start_quality,
                     FastqType ftype,
                     int depth)
{
    fprintf(out_fh, "%s\t%i\t%c\t%i\t",
            reference, position, reference_base, depth);
    
    size_t ftype_offset = FastqTypeOffset(ftype);

    fwrite(start_base, sizeof(start_base[0]), depth, out_fh);
    fprintf(out_fh, "\t");
    for (int si = 0; si != depth; ++si)
    {
        fputc(QualityToQualityCode(start_quality[si], ftype_offset), out_fh);
    }
    fprintf(out_fh, "\n");
}



//read in an unknown length line until a newline
char * read_unknown_length_line(FILE *file, int initial_buffer_size) 
{
    int current_size = initial_buffer_size;
    int char_num = 0; 
    char * line_buffer = 
        static_cast<char *>(malloc(current_size * sizeof(char)));

    int ch;
    for (ch=fgetc(file); ch != EOF && ch != '\n'; ch=fgetc(file)) 
    {
        if (char_num + 2 > current_size) 
        {
            current_size = 2 * char_num + 2; 
            line_buffer = static_cast<char *>(realloc(line_buffer,current_size));
        }
        line_buffer[char_num++] = ch; 
    }

    line_buffer[char_num] = 0;
    return line_buffer;
}


//parse a file with space-separated numbers, allocating a pointer
//of doubles and returning the number of numbers parsed
double * ParseNumbersFile(char const* numbers_file, size_t * num_numbers)
{
    FILE * numbers_fh = fopen(numbers_file, "r");
    if (numbers_fh == NULL)
    {
        fprintf(stderr, "ParseNumbersFile: Couldn't open numbers file %s\n", 
                numbers_file);
        exit(1);
    }

    size_t initial_size = 256;
    size_t current_size = initial_size;

    double * numbers = new double[current_size];

    size_t b = 0;
    while (! feof(numbers_fh))
    {
        if (b + 2 > current_size)
        {
            current_size = 2 * b + 2; 
            numbers = static_cast<double *>
                (realloc(numbers, sizeof(double) * current_size));
        }
        fscanf(numbers_fh, "%lf ", &numbers[b]);
        ++b;
    }
    fclose(numbers_fh);
    *num_numbers = b;
    return numbers;

}





// Given a range [begin, end) and a number of chunks num_chunks to
// divide into, return the begin or end of chunk chunk, in [0,
// num_chunks), whether to return the begin or end of the range is
// determined by 'give_begin'
size_t range_chunk_offset(size_t begin, size_t end, size_t num_chunks, size_t chunk, bool give_begin)
{
    assert(begin < end);
    assert(chunk < num_chunks);
    return (end - begin) * (chunk + (give_begin ? 0 : 1)) / num_chunks;
}
