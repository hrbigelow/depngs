#ifndef _TOOLS_H
#define _TOOLS_H

#include <cstdlib>
#include <cstdio>
#include <cassert>

#include <vector>
#include <cfloat>
#include <cmath>
#include <map>


template <typename ReturnVal, typename Iterator>
    ReturnVal GetMapData(Iterator const& iterator)
{
    return (*iterator).second;
}


template <typename Iterator>
typename Iterator::value_type Dereference(Iterator const& iterator)
{
    return *iterator;
}


//Create a lookup table to go from char to index, given
//an index-to-char map expressed in a char *
void CharLookupTable(char const* index_to_char_map, int default_value, 
                     int * lookup_table);


//compute log(n1 + n2), where n1 or n2 are underflow numbers expressed
//in log base 2.  If the difference between log_n1 and log_n2 is
//greater than 32, the lesser number is treated as -infinity
float SumLog2NumbersTruncated(float log_n1, float log_n2);


//computes the log(n1 + ... + nN) where we're given an array of numbers
//expressed as {2^n1, 2^n2, ... , 2^nN}
template <typename Iterator>
float SumLog2Truncated(Iterator begin,
                       Iterator end,
                       float (*unpack)(Iterator const&))
{

    //scan for max value
    float max_val = -FLT_MAX; //initialize to lowest value at start
    for (Iterator n_iter = begin; n_iter != end; ++n_iter)
    {
        max_val = std::max(unpack(n_iter), max_val);
    }
    
    //adjust all values by a constant factor (up or down) so the max is zero
    size_t size = std::distance(begin, end);
    std::vector<float> log_numbers_adjust(size);
    float sum_numbers = 0.0;

    std::vector<float>::iterator log_number_adjust = log_numbers_adjust.begin();

    for (Iterator log_number = begin; log_number != end; ++log_number)
    {
        *log_number_adjust = unpack(log_number) - max_val;
        sum_numbers += std::max(exp2f(*log_number_adjust), 0.0f);
        ++log_number_adjust;
    }
    
    return log2f(sum_numbers) + max_val;

}


//calculate the running sum of numbers
//log_numbers = { log(1), log(2), log(5), log(4) }
//calculate     { log(1), log(1+2), log(1+2+5), log(1+2+5+4) }
//print out a cumulative histogram (log space) for a sampled posterior
//estimate of the fraction of true majority base in the sample given
//the observed bases and quals


//creates an output container of elements that are the result of
//one-by-one conversion of input elements, using map_function
//to convert them
template <typename InputContainer, typename OutputContainer>
    OutputContainer Map(typename InputContainer::iterator begin,
                        typename InputContainer::iterator end,
                        typename OutputContainer::value_type 
                        (*map_function)(typename InputContainer::value_type const&))
{
    OutputContainer output;
    typename OutputContainer::iterator out_start = output.begin();

    for (typename InputContainer::iterator input = begin; input != end; ++input)
    {
        output.insert(output.end(), map_function(*input));
    }
    return output;
}


//creates an output container of elements that are the result of
//one-by-one conversion of input elements, using a map_fold_function
//the map_fold_function accepts an 'accumulator' and an input element
//and returns a pair of new accumulator and output element.
template <typename InputContainer, typename OutputContainer, typename Accumulator>
    OutputContainer MapFold(Accumulator const& initial_accumulator,
                            typename InputContainer::iterator begin,
                            typename InputContainer::iterator end,
                            std::pair<Accumulator, typename OutputContainer::value_type>
                            (*map_function)(Accumulator const&, 
                                            typename InputContainer::value_type const&))
{
    OutputContainer output;
    typename OutputContainer::iterator out_start = output.begin();
    
    Accumulator accumulator = initial_accumulator;

    for (typename InputContainer::iterator input = begin;
         input != end; ++input)
    {
        std::pair<Accumulator, typename OutputContainer::value_type>
            next_value = map_function(accumulator, *input);
        accumulator = next_value.first;
        output.insert(output.end(), next_value.second);
    }
    return output;
}


//creates an output container of elements that are the result of
//conversion of each consecutive pair of input elements, using a map_fold_function
//the map_fold_function accepts an 'accumulator' and two input elements
//and returns a pair of new accumulator and output element.
template <typename InputContainer, typename OutputContainer, typename Accumulator>
    OutputContainer MapFold2(Accumulator const& initial_accumulator,
                             typename InputContainer::iterator current,
                             typename InputContainer::iterator const& end,
                             std::pair<Accumulator, typename OutputContainer::value_type>
                             (*map_function)(Accumulator const&, 
                                             typename InputContainer::value_type const&,
                                             typename InputContainer::value_type const&))
{
    OutputContainer output;
    if (std::distance(current, end) < 2)
    {
        return output;
    }

    typename OutputContainer::iterator out_start = output.begin();
    
    Accumulator accumulator = initial_accumulator;

    typename InputContainer::iterator prev = current;
    for (++current; current != end; ++current)
    {
        std::pair<Accumulator, typename OutputContainer::value_type>
            next_value = map_function(accumulator, *prev, *current);
        accumulator = next_value.first;
        output.insert(output.end(), next_value.second);
        prev = current;
    }
    return output;
}


//a map of sample points (x,y) along a function

//a map of sample points (x,y) along a function

typedef std::map<double, double> SAMPLE;

//a single sample point from a SAMPLE (see above)
typedef SAMPLE::value_type SAMPLE_POINT;


SAMPLE Log2Accumulate(SAMPLE::iterator begin, SAMPLE::iterator end);

//produce a new set of sample points at the same X values, representing
//an upper bound on the Norm
SAMPLE UpperBoundNormLog2(SAMPLE::iterator begin, SAMPLE::iterator end);


//produce a new set of sample points at the same X values, representing
//the cumulative area of the lower bound step function
SAMPLE LowerBoundNormLog2(SAMPLE::iterator begin, SAMPLE::iterator end);


char QualityToQualityCode(int quality, size_t quality_code_offset);

size_t QualityCodeToQuality(char quality_code, size_t quality_code_offset);


enum FastqType {
    Sanger, // 33
    Solexa, // 64
    Illumina13, // 64
    Illumina15, // 64
    Illumina18, // 33
    None // 0
};


size_t FastqTypeOffset(FastqType ftype);


//determine the fastq type from the string of quality_codes
FastqType get_fastq_type(char const* quality_codes);


float QualityToErrorProb(int quality);


int error_prob_to_quality(float error_prob);

void PrintPileupLine(FILE * out_fh,
                     char const* reference,
                     int position,
                     char reference_base,
                     char const* start_base,
                     int const* start_quality,
                     FastqType ftype,
                     int depth);

//read in an unknown length line until a newline, returning a freshly allocated
//buffer
char * read_unknown_length_line(FILE *file, int initial_buffer_size);






//find the x value for which curve(x) is within tolerance of target_value
//assume that curve is monotonic increasing from start to end
//here, start may be greater than end.
template <typename Curve>
float BinarySearchFunctional(float start, float end, float target_value,
                             float x_tolerance, float y_tolerance, Curve & curve)
{
    
    assert(x_tolerance > 0.0 && y_tolerance > 0.0);

    float middle_x = start;
    float middle_y = curve(middle_x);
    float half = end - start;

    while (fabsf(half) > x_tolerance && fabsf(middle_y - target_value) > y_tolerance)
    {
        half /= 2.0;
        middle_x = start + half;
        middle_y = curve(middle_x);
        if (middle_y < target_value)
        {
            start = middle_x;
        }
    }

    return middle_x;

}


//parse a file with space-separated numbers, allocating a pointer
//of doubles and returning the number of numbers parsed
double * ParseNumbersFile(char const* numbers_file, size_t * num_numbers);


// Given a range [begin, end) and a number of chunks num_chunks to
// divide into, return the begin or end of chunk chunk, in [0,
// num_chunks), whether to return the begin or end of the range is
// determined by 'give_begin'
size_t range_chunk_offset(size_t begin, size_t end, size_t num_chunks, size_t chunk, bool give_begin);


#endif //_TOOLS_H
