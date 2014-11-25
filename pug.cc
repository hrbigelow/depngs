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
            "-b INT      size of output buffer [8e6]\n"
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
#define MAX(a,b) ((a) < (b) ? (b) : (a))

/* taken from git source code */
#define alloc_nr(x) (((x)+16)*3/2)

/*
 * Realloc the buffer pointed at by variable 'x' so that it can hold
 * at least 'nr' entries; the number of entries currently allocated
 * is 'alloc', using the standard growing factor alloc_nr() macro.
 *
 * DO NOT USE any expression with side-effect for 'x', 'nr', or 'alloc'.
 */
#define ALLOC_GROW(x, nr, alloc)                    \
    do {                                            \
        if ((nr) > alloc) {                         \
            if (alloc_nr(alloc) < (nr))             \
                alloc = (nr);                       \
            else                                    \
                alloc = alloc_nr(alloc);            \
            x = realloc((x), alloc * sizeof(*(x))); \
        }                                           \
    } while (0)


struct query_range {
    unsigned beg_contig, beg_pos, end_contig, end_pos;
};


int less_query_range(const void *pa, const void *pb)
{
    const struct query_range
        *a = (struct query_range *)pa,
        *b = (struct query_range *)pb;

    return 
        a-> < b ? -1 : a > b ? 1 
        : a->beg < b->beg ? -1 : a->beg > b->beg ? 1
        : a->end < b->end ? -1 : a->end > b->end ? 1
        : 0;
}



int main_pug(int argc, char ** argv)
{
    char c;
    size_t max_mem = 1024l * 1024l * 1024l;
    size_t max_pileup_line_size = 1e6;
    size_t outbuf_size = 8e6;

    while ((c = getopt(argc, argv, "m:l:b:")) >= 0)
    {
        switch(c)
        {
        case 'm': max_mem = static_cast<size_t>(atof(optarg)); break;
        case 'l': max_pileup_line_size = static_cast<size_t>(atof(optarg)); break;
        case 'b': outbuf_size = static_cast<size_t>(atof(optarg)); break;
        default: return pug_usage(); break;
        }
    }
    if (argc - optind != 3)
    {
        return pug_usage();
    }

    const char *pileup_file = argv[optind];
    const char *locus_file = argv[optind + 1];
    const char *contig_order_file = argv[optind + 2];

    std::map<char const*, size_t, ltstr> contig_order;
    std::map<char const*, size_t, ltstr>::iterator contig_iter;

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
    std::vector<char *> query_lines_tmp;
    char ** query;
    
    FILE *locus_fh = open_if_present(locus_file, "r");
    char contig[200];
    unsigned num_queries = 0, num_alloc = 10;
    struct query_range *queries = 
        (struct query_range *)malloc(num_alloc * sizeof(struct query_range)),
        *q;
    
    while (fscanf(locus_fh, "%s\t%u\t%u\n", contig, &q->beg, &q->end) == 3)
    {
        contig_iter = contig_order.find(contig);
        if (contig_iter == contig_order.end)
        {
            fprintf(stderr, "Error: Contig %s not listed in contig_order.rdb file\n", contig);
            exit(1);
        }
        q->contig_index = (*contig_iter).second;
        ++q;
        ++num_queries;
        ALLOC_GROW(q, num_queries, num_records);
    }   
    fclose(locus_fh);

    less_locus_position less_locus;
    equal_locus_position equal_locus;
    less_locus.contig_order = &contig_order;
    equal_locus.contig_order = &contig_order;

    qsort(queries, queries + num_queries, less_query_range);
    
    FILE *pileup_fh = open_if_present(pileup_file, "r");

    struct stat pileup_stat;
    fstat(fileno(pileup_fh), &pileup_stat);
    size_t total_pileup_size = pileup_stat.st_size;

    std::vector<char *>::iterator titer;
    char * const* qbeg;
    char * const* qend;
    qbeg = query;
   
    char *target_line = new char[max_pileup_line_size + 1];

    char *out_buf = new char[outbuf_size + 1];
    char *out_ptr = out_buf;
    char *out_end = out_buf + outbuf_size;

    // create a sparse index of offsets
    struct line_index
    {
        char line[50];
        size_t file_offset;
    };

    // one on the end for the end offset
    size_t index_chunk_size = 1e8;
    size_t index_size = MAX((total_pileup_size / index_chunk_size),1) + 1;
    fprintf(stderr, "Building index with %Zu entries.\n", index_size);
    fflush(stderr);
    line_index *index = new line_index[index_size];
    for (size_t i = 0; i != index_size - 1; ++i)
    {
        // throw away partial line
        fgets(target_line, max_pileup_line_size, pileup_fh);
        index[i].file_offset = ftell(pileup_fh);
        fscanf(pileup_fh, "%50c", index[i].line);
        fseek(pileup_fh, index_chunk_size * i, SEEK_SET);
        if (i != 0 && i % (index_size / 10) == 0)
        {
            fprintf(stderr, "Finished %Zu entries.\n", i);
            fflush(stderr);
        }
    }
    fprintf(stderr, "Finished building index.\n");
    fflush(stderr);

    // fill in the sentinel with a dummy value that is greater than the last
    char dummy_contig[50];
    size_t dummy_pos;
    sscanf(query[num_query - 1], "%s\t%zu", dummy_contig, & dummy_pos);
    sprintf(index[index_size - 1].line, "%s\t%Zu", dummy_contig, dummy_pos + 1);
    index[index_size - 1].file_offset = total_pileup_size;

    fseek(pileup_fh, 0, SEEK_SET);

    // find the next range of targets that will both fit in chunk_size and doesn't 
    // contain a gap greater than index_chunk_size
    size_t chunk_size = max_mem;
    size_t max_index_chunks = chunk_size / index_chunk_size;
    char *last_fragment;

    qbeg = query;
    size_t bi = 0, ei;


    // throw out all queries that occur before the beginning of the
    // index.
    while (*qbeg != NULL && less_locus(*qbeg, index[0].line))
    {
        ++qbeg;
    }

    
    while (*qbeg != NULL)
    {
        // find qend such that [qbeg, qend) fits within
        // max_index_chunks, and there is at least one target to
        // retrieve for each index chunk in the range

        // advance bi to the lower_bound position of qbeg
        while (less_locus(index[bi].line, *qbeg))
            ++bi;

        --bi;

        ei = bi + 1; // start with the minimal range that contains [qbeg, qend)
        qend = qbeg + 1;

        while (1)
        {
            // phase 1: extend qend until it doesn't fit within the current [bi, ei)
            while (*qend != NULL && less_locus(*qend, index[ei].line))
                ++qend;

            if (ei - bi < max_index_chunks)
            {
                // it is okay to extend number of index chunks

                // phase 2: if qend could fit within [bi, ei), extend ei.
                if (*qend != NULL && ei != index_size - 1 && less_locus(*qend, index[ei + 1].line))
                    ++ei;

                else break;
            }
            else break;
        }

        fprintf(stderr, "Using index range [%Zu to %Zu) of %Zu for %Zu query lines\n",
                bi, ei, index_size, std::distance(qbeg, qend));

        // [qbeg, qend) now should fit within [index[bi], index[ei])
        // note: this does NOT mean that qend < index[ei], just that qend - 1 < index[ei]

        // now, parse the whole chunk
        size_t this_chunk_size = index[ei].file_offset - index[bi].file_offset;
        char *chunk_buffer = new char[this_chunk_size + 1];
        
        fseek(pileup_fh, index[bi].file_offset, SEEK_SET);
        fread(chunk_buffer, 1, this_chunk_size, pileup_fh);
        std::vector<char *> target_lines = 
            FileUtils::find_complete_lines_nullify(chunk_buffer, & last_fragment);

        std::vector<char *>::iterator titer = target_lines.begin();
        assert(qbeg < qend);

        // print out any loci in [qbeg, qend) that exist in target_lines
        // we can break out of this loop either if all query loci are printed
        // or if we run out of target lines (which could happen if 
        while (qbeg != qend && titer != target_lines.end())
        {
            // advance titer
            titer = std::lower_bound(titer, target_lines.end(), *qbeg, less_locus);

            // advance qbeg
            while (qbeg != qend && less_locus(*qbeg, *titer))
            {
                ++qbeg;
            }

            if (qbeg == qend)
            {
                break;
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
                ++qbeg;
            }

        }
        qbeg = qend;

        delete chunk_buffer;

    }

    // write out remaining cached lines
    write(1, out_buf, out_ptr - out_buf);
    fsync(1);
    
    for (std::map<char const*, size_t, ltstr>::iterator cit = contig_order.begin();
         cit != contig_order.end(); ++cit)
    {
        delete (*cit).first;
    }

    delete locus_buf;
    delete target_line;
    delete out_buf;
    delete index;

    free(queries);

    fclose(pileup_fh);

    return 0;
}
