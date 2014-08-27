/*
From: http://www.1000genomes.org/wiki/Analysis/Variant%20Call%20Format/vcf-variant-call-format-version-41

##fileformat=VCFv4.1
##fileDate=20090805
##source=myImputationProgramV3.1
##reference=file:///seq/references/1000GenomesPilot-NCBI36.fasta
##contig=<ID=20,length=62435964,assembly=B36,md5=f126cdf8a6e0c7f379d618ff66beb2da,species="Homo sapiens",taxonomy=x>
##phasing=partial
##INFO=<ID=NS,Number=1,Type=Integer,Description="Number of Samples With Data">
##INFO=<ID=DP,Number=1,Type=Integer,Description="Total Depth">
##INFO=<ID=AF,Number=A,Type=Float,Description="Allele Frequency">
##INFO=<ID=AA,Number=1,Type=String,Description="Ancestral Allele">
##INFO=<ID=DB,Number=0,Type=Flag,Description="dbSNP membership, build 129">
##INFO=<ID=H2,Number=0,Type=Flag,Description="HapMap2 membership">
##FILTER=<ID=q10,Description="Quality below 10">
##FILTER=<ID=s50,Description="Less than 50% of samples have data">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype Quality">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read Depth">
##FORMAT=<ID=HQ,Number=2,Type=Integer,Description="Haplotype Quality">
#CHROM POS     ID        REF    ALT     QUAL FILTER INFO                              FORMAT      NA00001        NA00002        NA00003
20     14370   rs6054257 G      A       29   PASS   NS=3;DP=14;AF=0.5;DB;H2           GT:GQ:DP:HQ 0|0:48:1:51,51 1|0:48:8:51,51 1/1:43:5:.,.
20     17330   .         T      A       3    q10    NS=3;DP=11;AF=0.017               GT:GQ:DP:HQ 0|0:49:3:58,50 0|1:3:5:65,3   0/0:41:3
20     1110696 rs6040355 A      G,T     67   PASS   NS=2;DP=10;AF=0.333,0.667;AA=T;DB GT:GQ:DP:HQ 1|2:21:6:23,27 2|1:2:0:18,2   2/2:35:4
20     1230237 .         T      .       47   PASS   NS=3;DP=13;AA=T                   GT:GQ:DP:HQ 0|0:54:7:56,60 0|0:48:4:51,51 0/0:61:2
20     1234567 microsat1 GTC    G,GTCT  50   PASS   NS=3;DP=9;AA=G                    GT:GQ:DP    0/1:35:4       0/2:17:2       1/1:40:3


CHROM  -- as given
POS    -- as given
ID     -- leave blank
REF    -- as given
ALT    -- only handle single-base substitutions.
QUAL   -- not sure how there can be only a single quality here, given that the line pertains to multple samples
FILTER -- in this case, just 'PASS', because we want to annotate all of these.
INFO   -- if possible, just leave this field blank.  Hopefully, this is still a valid VCF file
FORMAT -- GT is the only required field, and it will be unphased.  (See below for other possible fields to include)

Other format fields:

GQ:  Since we are only using this as input to an annotator, this should be unnecessary
DP:  Ditto
HQ:  Ditto -- in fact, don't know what this is...


For genotype, if the sample has sample points, use the mean of sample
points as the best estimate for composition.  Then, take the closest
diploid composition as the diploid representation.

After compiling all diploid representations, determine the presence of
any ALT alleles. (this should be easy).

Then, output the Format field for each of these.

Questions about the VCF format in general:

How is the 'QUAL' field derived?  What is the formula?

INFO fields

* How is the AA (Ancestral Allele) value determined?
* IS the NS (Number of samples with data) field necessary?
* What is the DP (Total Depth) field used for?
* What is the AF (Allele Frequency) field used for?

Look at the sample_details field to get the composition estimate.
look at is_next to tell whether this is the next one or not.  If not,
it means there is no data for this locus for this sample, and you then
need to just output an arbitrary genotype.


 */

#include <gsl/gsl_math.h>

#include "vcf.h"
#include "dist_worker.h"
#include "pileup_tools.h"
#include "comp_functor.h"

double diploid_points[10][4] = {
    { 1.0, 0.0, 0.0, 0.0 },
    { 0.0, 1.0, 0.0, 0.0 },
    { 0.0, 0.0, 1.0, 0.0 },
    { 0.0, 0.0, 0.0, 1.0 },
    { 0.5, 0.5, 0.0, 0.0 },
    { 0.5, 0.0, 0.5, 0.0 },
    { 0.5, 0.0, 0.0, 0.5 },
    { 0.0, 0.5, 0.5, 0.0 },
    { 0.0, 0.5, 0.0, 0.5 },
    { 0.0, 0.0, 0.5, 0.5 }
};

#define MIN(a,b) ((a) < (b) ? (a) : (b))
#define MAX(a,b) ((a) < (b) ? (b) : (a))


// return an upper bound on the length of a vcf line with num_samples
size_t vcf_locus_bytes(size_t num_samples)
{
    return 10 + 16 + 1 + 1 + 5 + 2 + 4 + 0 + 2 + (num_samples * 3)
        + (9 + num_samples);
}


/*
  write the next vcf line to outbuf.  return next write position.
  
 */
char * print_vcf_line(sample_details * samples,
                      dist_worker_input * input,
                      char * outbuf)
{
    // find the first sample for which is_next is true
    PileupSummary * locus = NULL;
    for (size_t s = 0; s != input->num_samples; ++s)
    {
        if (samples[s].is_next)
        {
            locus = samples[s].locus;
            break;
        }
    }
    if (locus == NULL)
    {
        fprintf(stderr, "Error: print_vcf_line: no sample has the is_next flag set\n");
        exit(10);
    }

    outbuf += sprintf(outbuf, 
                      "%s\t%i\t.\t%c",
                      locus->reference,
                      locus->position,
                      locus->reference_base);

    bool alts[4] = { false, false, false, false };
    int * min_inds = new int[input->num_samples];

    for (size_t s = 0; s != input->num_samples; ++s)
    {
        /*
          1. get (A,C,G,T) composition estimate from samples[s], or NULL estimate if ! is_next
          2. compute nearest diploid genotype (among the 10)
          3. set a flag to indicate the presence of one of the alts
         */
        if (samples[s].is_next)
        {
            double mean_point[4] = { 0.0, 0.0, 0.0, 0.0 };
            const size_t np = samples[s].num_sample_points;

            // compute mean
            double * point = samples[s].sample_points;
            for (size_t p = 0; p != np; ++p)
            {
                mean_point[0] += point[0];
                mean_point[1] += point[1];
                mean_point[2] += point[2];
                mean_point[3] += point[3];
                point += 4;
            }
            mean_point[0] /= static_cast<double>(np);
            mean_point[1] /= static_cast<double>(np);
            mean_point[2] /= static_cast<double>(np);
            mean_point[3] /= static_cast<double>(np);

            // find nearest diploid point
            double min_dist = 3.0; // (this is larger than the maximum possible distance in base space
            min_inds[s] = 0;
            for (size_t p = 0; p != 10; ++p)
            {
                double dist2 =
                    gsl_pow_2(mean_point[0] * diploid_points[p][0])
                    + gsl_pow_2(mean_point[1] * diploid_points[p][1])
                    + gsl_pow_2(mean_point[2] * diploid_points[p][2])
                    + gsl_pow_2(mean_point[3] * diploid_points[p][3]);
                if (dist2 < min_dist)
                {
                    min_dist = dist2;
                    min_inds[s] = p;
                }
            }
            
            alts[0] = alts[0] || (diploid_points[min_inds[s]][0] > 0);
            alts[1] = alts[1] || (diploid_points[min_inds[s]][1] > 0);
            alts[2] = alts[2] || (diploid_points[min_inds[s]][2] > 0);
            alts[3] = alts[3] || (diploid_points[min_inds[s]][3] > 0);

        }
        else
        {
            // indicates to use a default dummy value later.
            min_inds[s] = -1;
        }
    }
    
    // now, min_inds and alts are initialized
    char bases[] = "ACGTN";
    int allele_ind[4];
    size_t ref_ind = strchr(bases, toupper(locus->reference_base)) - bases;
    if (ref_ind < 4)
    {
        allele_ind[ref_ind] = 0;
        alts[ref_ind] = false;
    }

    int allele_num = 1;
    char alt_list[20];
    alt_list[0] = '\0';
    char sep[2];
    sep[0] = '\0';
    sep[1] = '\0';

    for (size_t ind = 0; ind != 4; ++ind)
    {
        if (alts[ind])
        {
            allele_ind[ind] = allele_num++;
            strncat(alt_list, sep, 1);
            strncat(alt_list, bases + ind, 1);
            sep[0] = ',';
        }
    }

    size_t vcf_qual = 0; // dummy value
    char vcf_info[] = "";
    char vcf_format[] = "GT";

    outbuf += sprintf(outbuf,
                      "\t%s\t%Zu\tPASS\t%s\t%s",
                      alt_list,
                      vcf_qual,
                      vcf_info,
                      vcf_format);

    for (size_t s = 0; s != input->num_samples; ++s)
    {
        int locus_allele[2] = { 0, 0 };
        if (min_inds[s] != -1)
        {
            int a = 0;
            for (size_t b = 0; b != 4; ++b)
            {
                double p = diploid_points[min_inds[s]][b];
                while (p > 0)
                {
                    locus_allele[a++] = allele_ind[b];
                    p -= 0.5;
                }
            }
        }
        outbuf += sprintf(outbuf, 
                          "\t%i/%i", 
                          MIN(locus_allele[0], locus_allele[1]), 
                          MAX(locus_allele[0], locus_allele[1]));
    }

    outbuf += sprintf(outbuf, "\n");

    delete min_inds;

    return outbuf;
}
