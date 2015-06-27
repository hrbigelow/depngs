/* Convert one BAM file to a pileup file, using a streaming,
   multi-threaded approach. */

int pileup_usage()
{
    char *tmp_require = bam_flag2str(mplp->rflag_require);
    char *tmp_filter  = bam_flag2str(mplp->rflag_filter);
    
    // Display usage information, formatted for the standard 80 columns.
    // (The unusual string formatting here aids the readability of this
    // source code in 80 columns, to the extent that's possible.)

    fprintf(stderr,
            "\nUsage: dep pileup [options] in.bam out.pileup\n\n"
            "Options:\n"
            "-t  INT    number of threads [1]\n"
            "-F  INT    Output Phred Quality Encoding offset (33 or 64) [33]\n"
            "-A  FLAG   do not discard anomalous read pairs\n"
            "-B  FLAG   disable BAQ (per-Base Alignment Quality)\n"
            "-C  INT    adjust mapping quality; recommended:50, disable:0 [0]\n"
            "-E  FLAG   recalculate BAQ on the fly, ignore existing BQs\n"
            "-f  STR    faidx indexed reference sequence filename\n"
            "-G  STR    file list of read groups to exclude\n"
            "-l  STR    file with ranges (chr start end) to process\n"
            "-q  INT    min mapQ to include alignments [%d]\n"
            "-Q  INT    min baseQ/BAQ to include in pileup [%d]\n"
            "-x  STR    flags that must be present for inclusion (string or int) [%s]\n"
            "-X  STR    flags that must be absent for inclusion (string or int) [%s]\n",
            mplp->min_mq,
            mplp->min_baseQ,
            tmp_require,
            tmp_filter);
    
    free(tmp_require);
    free(tmp_filter);
}

int main_pileup(int argc, char **argv)
{
    char c;
    
    /* adapted from samtools/bam_plcmd.c.*/
    while ((c = getopt(argc, argv, "t:F:ABC:Ef:G:l:q:Q:x:X:")) >= 0) {
        switch (c) {
        case 't': n_threads = strtol_errmsg(optarg, "-t (n_threads)"); break;
        case 'F': phred_offset = strtol_errmsg(optarg, "-F (phred_offset)"); break;
        case 'A': use_orphan = 1; break;
        case 'B': mplp.flag &= ~MPLP_REALN; break;
        case 'C': mplp.capQ_thres = atoi(optarg); break;
        case 'E': mplp.flag |= MPLP_REDO_BAQ; break;
        case 'f':
            mplp.fai = fai_load(optarg);
            if (mplp.fai == 0) return 1;
            mplp.fai_fname = optarg;
            break;
        case 'G': {
                FILE *fp_rg;
                char buf[1024];
                mplp.rghash = khash_str2int_init();
                if ((fp_rg = fopen(optarg, "r")) == 0)
                    fprintf(stderr, "(%s) Fail to open file %s. Continue anyway.\n", __func__, optarg);
                while (!feof(fp_rg) && fscanf(fp_rg, "%s", buf) > 0) // this is not a good style, but forgive me...
                    khash_str2int_inc(mplp.rghash, strdup(buf));
                fclose(fp_rg);
            }
            break;
        case 'l':
                  mplp.bed = bed_read(optarg);
                  if (!mplp.bed) { print_error_errno("Could not read file \"%s\"", optarg); return 1; }
                  break;
        case 'q': mplp.min_mq = atoi(optarg); break;
        case 'Q': mplp.min_baseQ = atoi(optarg); break;

        case 'x' :
            mplp.rflag_require = bam_str2flag(optarg);
            if ( mplp.rflag_require<0 ) { fprintf(stderr,"Could not parse -x %s\n", optarg); return 1; }
            break;
        case 'X' :
            mplp.rflag_filter = bam_str2flag(optarg);
            if ( mplp.rflag_filter<0 ) { fprintf(stderr,"Could not parse --X %s\n", optarg); return 1; }
            break;

        case  4 : mplp.openQ = atoi(optarg); break;
        case 'p': mplp.flag |= MPLP_PER_SAMPLE; break;
        case 'D': mplp.fmt_flag |= B2B_FMT_DP; fprintf(stderr, "[warning] samtools mpileup option `-D` is functional, but deprecated. Please switch to `-t DP` in future.\n"); break;
        case 'S': mplp.fmt_flag |= B2B_FMT_SP; fprintf(stderr, "[warning] samtools mpileup option `-S` is functional, but deprecated. Please switch to `-t SP` in future.\n"); break;
        case 'V': mplp.fmt_flag |= B2B_FMT_DV; fprintf(stderr, "[warning] samtools mpileup option `-V` is functional, but deprecated. Please switch to `-t DV` in future.\n"); break;

        default:
            fprintf(stderr,"Invalid option: '%c'\n", c);
            return 1;
        }
    }

    if (argc - optind != 2) return pileup_usage();
        
    char *bam_file = argv[optind];
    char *plp_file = argv[optind + 1];
    /* */


}
