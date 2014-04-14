#include <cstdio>
#include <unistd.h>
#include <cstdlib>
#include <map>
#include <algorithm>

#include <sys/stat.h>

#include "locus_comp.h"
#include "tools.h"
#include "../samutil/file_utils.h"

int pug_usage()
{
    fprintf(stderr,
            "\nUsage: dep pug [options] sample.pileup loci_to_retrieve.rdb contig_order.rdb\n"
            "Options:\n\n"
            "-m INT      number bytes of memory to use [1e9]\n"
            "-l INT      maximum length of a pileup line in bytes.  {Nuisance parameter} [100000]\n"
            "\n"
            "loci_to_retrieve.rdb has lines like:\n"
            "chr1<tab>19583984\n"
            "chr1<tab>19598940\n"
            "...\n"
            "It need not be in order.  Also, it is okay if loci do not appear in the sample.pileup file\n"
            "\n"
            "contig_order.rdb has lines of <contig><tab><ordering>\n"
            "It must be consistent and complete with the orderings of the contigs mentioned in\n"
            "all pileup input files\n"
            );
    return 1;
}

#define MIN(a,b) ((a) < (b) ? (a) : (b))

int main_pug(int argc, char ** argv)
{
    char c;
    size_t max_mem = 1024l * 1024l * 1024l;
    size_t max_pileup_line_size = 1e6;

    while ((c = getopt(argc, argv, "m:l:")) >= 0)
    {
        switch(c)
        {
        case 'm': max_mem = static_cast<size_t>(atof(optarg)); break;
        case 'l': max_pileup_line_size = static_cast<size_t>(atof(optarg)); break;
        default: return pug_usage(); break;
        }
    }
    if (argc - optind != 3)
    {
        return pug_usage();
    }

    char const* pileup_file = argv[optind];
    char const* locus_file = argv[optind + 1];
    char const* contig_order_file = argv[optind + 2];

    std::map<char const*, size_t, ltstr> contig_order;
    // parse contig_order file
    {
        FILE * contig_order_fh = open_if_present(contig_order_file, "r");
        while (! feof(contig_order_fh))
        {
            char * contig = new char[100];
            size_t order;
            fscanf(contig_order_fh, "%s\t%zu\n", contig, & order);
            contig_order.insert(std::make_pair(contig, order));
        }
        fclose(contig_order_fh);
    }

    // parse locus file
    char * locus_buf = NULL;
    std::vector<char *> query_lines;
    {
        FILE * locus_fh = open_if_present(locus_file, "r");
        struct stat locus_stat;
        fstat(fileno(locus_fh), &locus_stat);
        size_t total_size = locus_stat.st_size;
        if (total_size > max_mem)
        {
            fprintf(stderr, "Error: total size of the loci_to_retrieve.rdb is %Zu,\n"
                    "larger than -m option (max memory) of %Zu.  It needs to be loaded\n"
                    "in memory however.  Please increase -m\n",
                    total_size, max_mem);
            exit(10);
        }
        locus_buf = new char[total_size];
        fread(locus_buf, 1, total_size, locus_fh);
        fclose(locus_fh);
        char * last_fragment;
        query_lines = FileUtils::find_complete_lines_nullify(locus_buf, &last_fragment);
        
    }

    less_locus_position less_locus;
    equal_locus_position equal_locus;
    less_locus.contig_order = &contig_order;
    equal_locus.contig_order = &contig_order;

    std::sort(query_lines.begin(), query_lines.end(), less_locus);
    
    FILE * pileup_fh = open_if_present(pileup_file, "r");

    struct stat pileup_stat;
    fstat(fileno(pileup_fh), &pileup_stat);
    size_t total_pileup_size = pileup_stat.st_size;

    std::vector<char *>::iterator qcur, qbeg, qend, titer;
    qbeg = query_lines.begin();
    // size_t min_loci_for_linear = 1;
    // size_t average_line_length = 100;
    // size_t jump = min_loci_for_linear * average_line_length;

    // char query_locus[100];
    // char query_contig[100];
    // size_t query_position;
    // size_t fread_elapsed_nsec;
    char * target_line = new char[max_pileup_line_size + 1];

    size_t outbuf_size = max_mem;

    char * out_buf = new char[outbuf_size + 1];
    char * out_ptr = out_buf;
    char * out_end = out_buf + outbuf_size;

    // size_t lb, ub, cur;

    // size_t fseeks = 0;
    // size_t fgetss = 0;
    // size_t freads = 0;

    // compute contig sizes
    

    // create a sparse index of offsets
    struct line_index
    {
        char line[50];
        size_t file_offset;
    };

    // one on the end for the end offset
    size_t index_chunk_size = 1e8;
    size_t index_size = (total_pileup_size / index_chunk_size) + 1;
    line_index * sparse_index = new line_index[index_size];
    for (size_t i = 0; i != index_size - 1; ++i)
    {
        // throw away partial line
        fgets(target_line, max_pileup_line_size, pileup_fh);
        sparse_index[i].file_offset = ftell(pileup_fh);
        fscanf(pileup_fh, "%50c", sparse_index[i].line);
        fseek(pileup_fh, index_chunk_size * i, SEEK_SET);
    }

    // fill in the sentinel with a dummy value that is greater than the last
    char dummy_contig[50];
    size_t dummy_pos;
    sscanf(*query_lines.rbegin(), "%s\t%zu", dummy_contig, & dummy_pos);
    sprintf(sparse_index[index_size - 1].line, "%s\t%Zu", dummy_contig, dummy_pos + 1);
    sparse_index[index_size - 1].file_offset = total_pileup_size;

    fseek(pileup_fh, 0, SEEK_SET);

    // find the next range of targets that will both fit in chunk_size and doesn't 
    // contain a gap greater than index_chunk_size
    size_t chunk_size = max_mem;
    size_t max_index_chunks = chunk_size / index_chunk_size;
    char * last_fragment;

    qbeg = query_lines.begin();
    size_t bi = 0, ei;
    
    while (qbeg != query_lines.end())
    {
        // find qend such that [qbeg, qend) fits within
        // max_index_chunks, and there is at least one target to
        // retrieve for each index chunk in the range

        // advance bi to the lower_bound position of qbeg
        while (less_locus(sparse_index[bi].line, *qbeg))
        {
            ++bi;
        }
        --bi;

        ei = bi + 1; // start with the minimal range that contains [qbeg, qend)
        qend = qbeg + 1;

        while (1)
        {
            // phase 1: extend qend until it doesn't fit within the current [bi, ei)
            while (qend != query_lines.end() && less_locus(*qend, sparse_index[ei].line))
            {
                ++qend;
            }

            if (ei - bi < max_index_chunks)
            {
                // it is okay to extend number of index chunks

                // phase 2: if qend could fit within [bi, ei), extend ei.
                if (qend != query_lines.end() && ei != index_size - 1 && less_locus(*qend, sparse_index[ei + 1].line))
                {
                    ++ei;
                }
                else
                {
                    break;
                }
            }
            else
            {
                break;
            }
        }

        fprintf(stderr, "Using index range [%Zu to %Zu) for %Zu query lines\n",
                bi, ei, std::distance(qbeg, qend));

        // now, parse the whole chunk
        size_t this_chunk_size = sparse_index[ei].file_offset - sparse_index[bi].file_offset;
        char * chunk_buffer = new char[this_chunk_size + 1];
        
        fseek(pileup_fh, sparse_index[bi].file_offset, SEEK_SET);
        fread(chunk_buffer, 1, this_chunk_size, pileup_fh);
        std::vector<char *> target_lines = FileUtils::find_complete_lines_nullify(chunk_buffer, & last_fragment);
        std::vector<char *>::iterator titer = target_lines.begin();
        while (qbeg != qend)
        {
            titer = std::lower_bound(titer, target_lines.end(), *qbeg, less_locus);

            while (qbeg != qend && less_locus(*qbeg, *titer))
            {
                ++qbeg;
            }
            
            if (equal_locus(*qbeg, *titer))
            {
                // print out this qbeg
                size_t ll = strlen(*titer);
                if (out_ptr + ll > out_end)
                {
                    // flush buffer
                    write(1, out_buf, out_ptr - out_buf);
                    fsync(1);
                    out_ptr = out_buf;
                }
                strcpy(out_ptr, *titer);
                out_ptr += ll;
                *out_ptr = '\n';
                ++out_ptr;
            }

            ++qbeg;
        }
        delete chunk_buffer;

    }

    // write out remaining cached lines
    write(1, out_buf, out_ptr - out_buf);
    fsync(1);
    
    // fprintf(stderr, "fseeks: %Zu, fgetss: %Zu, freads: %Zu\n",
    //         fseeks, fgetss, freads);

    for (std::map<char const*, size_t, ltstr>::iterator cit = contig_order.begin();
         cit != contig_order.end(); ++cit)
    {
        delete (*cit).first;
    }

    delete target_line;
    delete out_buf;
    fclose(pileup_fh);

    return 0;
}

    
    /*
    while (qbeg != query_lines.end())
    {
        // start binary seek-search.  at the end of this,
        // pileup_fh will be positioned at the beginning of the
        // first line that is not less than qbeg
        
        // we will either be at the beginning of a new line, or at the
        // end of a file
        lb = ftell(pileup_fh);
        ub = MIN(lb + jump, total_pileup_size);

        // phase 1: grow interval forward
        while (ub != total_pileup_size)
        {
            assert(lb < ub);
            fseek(pileup_fh, ub, SEEK_SET);
            fgets(target_line, max_pileup_line_size, pileup_fh);
            ++fseeks;
            ++fgetss;
            cur = ftell(pileup_fh);

            if (cur == total_pileup_size)
            {
                break;
            }
            else
            {
                fgets(target_line, max_pileup_line_size, pileup_fh);
                ++fgetss;
                if (less_locus(*qbeg, target_line))
                {
                    ub = cur;
                    break;
                }
                else
                {
                    size_t s = ub - lb;
                    lb = cur;
                    ub = MIN(cur + 2 * s, total_pileup_size);
                }
            }
        }

        // phase 2: shrink interval
        while (ub - lb > jump)
        {
            fseek(pileup_fh, (lb + ub) / 2, SEEK_SET);
            fgets(target_line, max_pileup_line_size, pileup_fh);
            ++fseeks;
            ++fgetss;
            cur = ftell(pileup_fh);

            if (cur == ub)
            {
                cur = lb;
                break;
            }
            else
            {
                fgets(target_line, max_pileup_line_size, pileup_fh);
                fseek(pileup_fh, cur, SEEK_SET);
                ++fseeks;
                ++fgetss;

                if (less_locus(*qbeg, target_line))
                {
                    ub = cur;
                }
                else if (equal_locus(*qbeg, target_line))
                {
                    break;
                }
                else
                {
                    lb = cur;
                }
            }
        }

        fseek(pileup_fh, lb, SEEK_SET);
        
        // 2. search for qend.  it is only used to grab a new chunk of the target file
        // sscanf((*qbeg), "%s\t%zu", query_contig, &query_position);
        // sprintf(query_locus, "%s\t%Zu", query_contig, query_position + min_loci_for_linear);
        // qend = std::upper_bound(qbeg, query_lines.end(), query_locus, less_locus);
        
        // 3. populate target_lines.
        size_t bytes_wanted = ub - lb;
        size_t max_bytes = bytes_wanted + max_pileup_line_size;
        char * target_buf = new char[max_bytes + 1];
        char * last_fragment;
        size_t bytes_read = 
            FileUtils::read_until_newline(target_buf, bytes_wanted, max_pileup_line_size,
                                          pileup_fh, & fread_elapsed_nsec);

        ++freads;
        std::vector<char *> target_lines = FileUtils::find_complete_lines_nullify(target_buf, & last_fragment);
        titer = target_lines.begin();

        // 4. traverse and output.  initially, titer is guaranteed to be not less than
        // qbeg.
        for (titer = target_lines.begin(); titer != target_lines.end(); ++titer)
        {

            while (qbeg != query_lines.end() && less_locus(*qbeg, *titer))
            {
                ++qbeg;
            }

            if (qbeg != query_lines.end() && equal_locus(*qbeg, *titer))
            {
                // target contains the locus that we want.  print it, optionally flushing the buffer first
                size_t ll = strlen(*titer);
                if (out_ptr + ll > out_end)
                {
                    // flush buffer
                    write(1, out_buf, out_ptr - out_buf);
                    fsync(1);
                    out_ptr = out_buf;
                }
                strcpy(out_ptr, *titer);
                out_ptr += ll;
                *out_ptr = '\n';
                ++out_ptr;
            }
            else
            {
                // target is missing the locus we want.  do nothing
            }
        }
        delete target_buf;
                                                              
    }
    */

