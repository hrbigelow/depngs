#ifndef _VCF_H
#define _VCF_H

#include "dist_worker.h"

size_t vcf_locus_bytes(size_t num_samples);

char * print_vcf_line(sample_details * samples, 
                      dist_worker_input * input,
                      char * outbuf);

#endif // _VCF_H
