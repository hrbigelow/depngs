//Convert fastq to fastq, with a different qual encoding
#include <algorithm>
#include <cstdio>

int main(int argc, char ** argv){

    if (argc != 3)
    {
        fprintf(stderr, 
                "Usage: fastq2fastq input.fastq output.fastq\n"
                "\n"
                "Converts to other offset fastq (there are only two unique ones)\n"
                "Assumes input fastq is one of these (from Wikipedia):\n"
                
                "\n"
                "  SSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSS.....................................................\n"
                "  ..........................XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX......................\n"
                "  ...............................IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII......................\n"
                "  .................................JJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJ......................\n"
                "  ..LLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLL....................................................\n"
                "  !\"#$%%&'()*+,-./0123456789:;<=>?@ABCDEFGHIJKLMNOPQRSTUVWXYZ[\\]^_`abcdefghijklmnopqrstuvwxyz{|}~\n"
                "  |                         |    |        |                              |                     |\n"
                " 33                        59   64       73                            104                   126\n"
                "\n"
                " S - Sanger        Phred+33,  raw reads typically (0, 40)\n"
                " X - Solexa        Solexa+64, raw reads typically (-5, 40)\n"
                " I - Illumina 1.3+ Phred+64,  raw reads typically (0, 40)\n"
                " J - Illumina 1.5+ Phred+64,  raw reads typically (3, 40)\n"
                " L - Illumina 1.8+ Phred+33,  raw reads typically (0, 41)\n"
                "\n");
        return 0;

    }

    char * input_file = argv[1];
    char * output_file = argv[2];

    FILE * input_fh = fopen(input_file, "r");
    FILE * output_fh = fopen(output_file, "w");

    if (input_fh == NULL)
    {
        fprintf(stderr, "Couldn't open input file %s\n", input_file);
        exit(1);
    }
    if (output_fh == NULL)
    {
        fprintf(stderr, "Couldn't open output file %s\n", output_file);
        exit(1);
    }

    //pre-check to see what kind of input file we have
    char * file_type_strings[] = { 
        "Sanger, Phred+33, range (0,40)",
        "Solexa, Solexa+64, range (-5,40)",
        "Illumina 1.3+, Phred+64, (0,40)",
        "Illumina 1.5+, Phred+64, (3,40)" 
    };

    int nfields_read;
    char id[1024];
    char sequence[1024];
    char spacer[1024];
    char quality_string[1024];
    
    char * qual;

    //try to eliminate other files by finding
    bool file_types[] = { true, true, true, true };
    char * min_qualcodes = "!;@B";
    char * bound_qualcodes = "Jiii";

    bool legal_values[4][256];
    for (size_t st = 0; st != 4; ++st)
    {
        for (size_t vt = 0; vt != 256; ++vt)
        {
            legal_values[st][vt] = 
                min_qualcodes[st] <= static_cast<char>(vt) 
                && static_cast<char>(vt) < bound_qualcodes[st];
        }
    }

    size_t file_type = 5;
    size_t chunks = 0;
    while (! feof(input_fh)){
        nfields_read = 
            fscanf(input_fh, "%*s\n%*s\n%*s\n%s\n", quality_string);
        
        ++chunks;
        for (qual = quality_string; *qual != '\0'; ++qual)
        {
            for (size_t qc = 0; qc != 4; ++qc)
            {
                file_types[qc] = file_types[qc] 
                    && legal_values[qc][static_cast<size_t>(*qual)];
            }
        }
        if (file_types[0] &&
            ! (file_types[1] || file_types[2] || file_types[3]))
        {
            file_type = 0;
            break;
        }
        else if ((file_types[1] || file_types[2] || file_types[3]) &&
                 ! file_types[0])
        {
            file_type = 1;
            break;
        }
        if (chunks > 1000000)
        {
            break;
        }
            
    }

    if (file_type == 5)
    {
        fprintf(stderr, "Couldn't determine input fastq file type\n");
        exit(1);
    }
    else
    {
        fprintf(stderr, "Input file is one of:\n");
        for (size_t f = 0; f != 4; ++f)
        {
            if (file_types[f])
            {
                fprintf(stderr, "%s\n", file_type_strings[f]);
            }
        }
        rewind(input_fh);
    }

    int input_output_diff = file_type == 0 ? (33 - 64) : (64 - 33);
    
    char qseq2fastq[256];
    for (int i=0; i < 256; ++i){
        qseq2fastq[i] = char(std::max(0,i - input_output_diff));
    }

    while (! feof(input_fh)){
        nfields_read = 
            fscanf(input_fh, "%s\n%s\n%s\n%s\n", id, sequence, spacer, quality_string);

        if (nfields_read != 4)
        {
            if (nfields_read > 0)
            {
                fprintf(stderr, "fastq2fastq: found badly formatted input group with %i lines.\n", nfields_read);
                exit(1);
            }
            break;
        }
        for (char * q = quality_string; *q != 0; ++q)
        { 
            *q = qseq2fastq[int(*q)];
        }

        fprintf(output_fh, "%s\n%s\n%s\n%s\n", id, sequence, spacer, quality_string);
    }
    fclose(input_fh);
    fclose(output_fh);
    return 0;
}
