#include <cstdlib>
#include <cstddef>
#include <cstdint>
#include <cstdio>
#include <algorithm>
#include <gsl/gsl_math.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_blas.h>
#include <pthread.h>
#include <cassert>
#include <vector>
#include <string.h>
#include <samutil/file_utils.h>

/*
  load an RDB file, and efficiently perform the pairwise query:

  select 
  a.snp_seq_changes_id,
  a.snp_sample_id sample1_id,
  b.snp_sample_id sample2_id, 
  sqrt(sum(power(a.comp - b.comp, 2))) dist_to_modes,
  sqrt(sum(power(case when a.
  from snp_mode_indiv a, snp_mode_indiv b, snp_pairwise_comp c
  where a.snp_seq_changes_id = b.snp_seq_changes_id
  and a.snp_sample_id = c.snp_sample1_id
  and b.snp_sample_id = c.snp_sample2_id
  group by a.snp_seq_changes_id, c.comparison_id
 */


/*
  input should be organized by (top to bottom), with cardinality given ()

  snp_sample_id (35)
  locus (11,000,000)
  hyp_base (4)
  stats (4) (rank, mode, lo, hi)

  output should be:

  comparison_id (N)
  locus (11,000,000)
  stats (3)

 */


/*
// find the nearest distance between two line segments in 3D
float nearest_dist(float * A1, float * A2, float * B1, float * B2)
{
    gsl_vector * vA1 = gsl_vector_alloc(3);
    gsl_vector * vA2 = gsl_vector_alloc(3);
    gsl_vector * vB1 = gsl_vector_alloc(3);
    gsl_vector * vB2 = gsl_vector_alloc(3);

    gsl_vector * vA = gsl_vector_alloc(3);
    gsl_vector * vB = gsl_vector_alloc(3);

    // 's' means 'same', 'x' means 'cross'
    // '1' means vector starts at A1, '2' means starts at A2
    // e.g. v1s means A1 to B1
    gsl_vector * v1s = gsl_vector_alloc(3);
    gsl_vector * v1x = gsl_vector_alloc(3);
    gsl_vector * v2s = gsl_vector_alloc(3);
    gsl_vector * v2x = gsl_vector_alloc(3);

    gsl_vector ** endpoints = { vA1, vA2, vB1, vB2 };
    gsl_vector ** struts = { v1s, v1x, v2s, v2x };
    float ** vals = { A1, A2, B1, B2 };
    float norm[4];
    double dot[8];

    // initialize the vectors
    for (size_t v = 0; v != 4; ++v)
    {
        for (size_t c = 0; c != 3; ++c)
        {
            gsl_vector_set(vecs[v], c, vals[v][c]);
        }
        norm[v] = gsl_blas_dnrm2(vecs[v]);
    }
    // initialize the line segment vectors
    gsl_vector_memcpy(vA, vA2);
    gsl_vector_memcpy(vB, vB2);

    gsl_vector_sub(vA, vA1);
    gsl_vector_sub(vB, vB1);

    // initialize the struts
    gsl_vector_memcpy(v1s, vB1);
    gsl_vector_memcpy(v1x, vB2);
    gsl_vector_memcpy(v2s, vB1);
    gsl_vector_memcpy(v2x, vB2);

    gsl_vector_sub(v1s, vA1);
    gsl_vector_sub(v1x, vA1);
    gsl_vector_sub(v2s, vA2);
    gsl_vector_sub(v2x, vA2);

     = gsl_vector_ddot(
    
}

float nearest_dist(float * seg1_start, float * seg1_end,
                   float * seg2_start, float * seg2_end)
{

    float const SMALL_NUM = 1e-30;
    float u[3];
    float v[3];
    float w[3];

    diff3(seg1_start, seg1_end, u);
    diff3(seg2_start, seg2_end, v);
    diff3(seg2_start, seg1_start, w);

    float a = dotp3(u,u);         // always >= 0
    float b = dotp3(u,v);
    float c = dotp3(v,v);         // always >= 0
    float d = dotp3(u,w);
    float e = dotp3(v,w);
    float D = a*c - b*b;        // always >= 0
    float sc, sN, sD = D;       // sc = sN / sD, default sD = D >= 0
    float tc, tN, tD = D;       // tc = tN / tD, default tD = D >= 0

    // compute the line parameters of the two closest points
    if (D < SMALL_NUM) { // the lines are almost parallel
        sN = 0.0;         // force using point P0 on segment S1
        sD = 1.0;         // to prevent possible division by 0.0 later
        tN = e;
        tD = c;
    }
    else {                 // get the closest points on the infinite lines
        sN = (b*e - c*d);
        tN = (a*e - b*d);
        if (sN < 0.0) {        // sc < 0 => the s=0 edge is visible
            sN = 0.0;
            tN = e;
            tD = c;
        }
        else if (sN > sD) {  // sc > 1  => the s=1 edge is visible
            sN = sD;
            tN = e + b;
            tD = c;
        }
    }

    if (tN < 0.0) {            // tc < 0 => the t=0 edge is visible
        tN = 0.0;
        // recompute sc for this edge
        if (-d < 0.0)
            sN = 0.0;
        else if (-d > a)
            sN = sD;
        else {
            sN = -d;
            sD = a;
        }
    }
    else if (tN > tD) {      // tc > 1  => the t=1 edge is visible
        tN = tD;
        // recompute sc for this edge
        if ((-d + b) < 0.0)
            sN = 0;
        else if ((-d + b) > a)
            sN = sD;
        else {
            sN = (-d +  b);
            sD = a;
        }
    }
    // finally do the division to get sc and tc
    sc = (fabsf(sN) < SMALL_NUM ? 0.0 : sN / sD);
    tc = (fabsf(tN) < SMALL_NUM ? 0.0 : tN / tD);

    // get the difference of the two closest points
    float dP[3];
    dP[0] = w[0] + (sc * u[0]) - (tc * v[0]);
    dP[1] = w[1] + (sc * u[1]) - (tc * v[1]);
    dP[2] = w[2] + (sc * u[2]) - (tc * v[2]);

    // norm(dP)
    return sqrt(dP[0] * dP[0] + dP[1] * dP[1] + dP[2] * dP[2]);
}
*/


struct input_stat
{
    uint32_t rank; // to ensure good alignment
    float mode;
    float lo;
    float hi;
    char base;
};


struct output_stat
{
    uint16_t sample1_id;
    uint16_t sample2_id;
    uint32_t locus_id;
    float dist_mode;
    float dist_lo;
    float dist_hi;
};

struct calc_input
{
    size_t locus_id;
    size_t sample1_id;
    size_t sample2_id;
    input_stat * sample1_beg;
    input_stat * sample1_end;
    input_stat * sample2_beg;
    output_stat * out_start;
};


struct comparison
{
    size_t sample1_id;
    size_t sample2_id;
};

inline float ssd(float a1, float a2, float a3, float a4,
                 float b1, float b2, float b3, float b4)
{
    return 
        gsl_pow_2(a1 - b1)
        + gsl_pow_2(a2 - b2)
        + gsl_pow_2(a3 - b3)
        + gsl_pow_2(a4 - b4);
}


/*
float dotp3(float * p1, float *p2)
{
    return p1[0] * p2[0]
        + p1[1] * p2[1]
        + p1[2] * p2[2];
}

void diff3(float * p1, float * p2, float * diff)
{
    diff[0] = p2[0] - p1[0];
    diff[1] = p2[1] - p1[1];
    diff[2] = p2[2] - p1[2];
}
*/


// compute p + s * g
void vecx(gsl_vector * p, gsl_vector * g, double s, gsl_vector * x)
{
    gsl_vector_memcpy(x, g);
    gsl_vector_scale(x, s);
    gsl_vector_add(x, p);
}


void disa(gsl_vector * x1, gsl_vector * x2, gsl_vector * da, double * lq)
{
    gsl_vector_memcpy(da, x2);
    gsl_vector_sub(da, x1);
    gsl_vector_mul(da, da);
    *lq = gsl_blas_dasum(da);
}


void difv(gsl_vector * p, gsl_vector * q, gsl_vector * r)
{
    gsl_vector_memcpy(r, p);
    gsl_vector_sub(r, q);
}


float nearest_dist(float * vp1, float * vq1, float * vp2, float * vq2)
{
    gsl_vector * p1, * q1, * p2, * q2;
    double a11, a12, a22, b1, b2, D0, D1, D2, s1, s2, n1, n2, sa, sb, dq, eq;
    gsl_vector * g1, * g2, * dp, * da;
    gsl_vector * x1, * x2;

    p1 = gsl_vector_alloc(3);
    q1 = gsl_vector_alloc(3);
    p2 = gsl_vector_alloc(3);
    q2 = gsl_vector_alloc(3);

    g1 = gsl_vector_alloc(3);
    g2 = gsl_vector_alloc(3);
    dp = gsl_vector_alloc(3);
    da = gsl_vector_alloc(3);

    x1 = gsl_vector_alloc(3);
    x2 = gsl_vector_alloc(3);

    gsl_vector * vecs[] = { p1, q1, p2, q2 };
    float * vals[] = { vp1, vq1, vp2, vq2 };

    for (size_t v = 0; v != 4; ++v)
    {
        for (size_t c = 0; c != 3; ++c)
        {
            gsl_vector_set(vecs[v], c, vals[v][c]);
        }
    }

    difv(q1, p1, g1);
    difv(q2, p2, g2);
    difv(p2, p1, dp);
    gsl_blas_ddot(g1, g1, & a11);
    gsl_blas_ddot(g1, g2, & a12); a12 *= -1.0;
    gsl_blas_ddot(g2, g2, & a22);
    gsl_blas_ddot(g1, dp, & b1);
    gsl_blas_ddot(g2, dp, & b2); b2 *= -1.0;

    // Cramer determinants
    D0 = a11 * a22 - a12 * a12;
    D1 = b1 * a22 - b2 * a12;
    D2 = -b1 * a12 + b2 * a11;

    // In-line distance
    int flag = 0;
    if (D0 > 0)
    {
        if (0 < D1 && D1 < D0) 
        {
            ++flag;
        }
        if (0 < D2 && D2 < D0)
        {
            ++flag;
        }
    }
    if (D0 < 0)
    {
        if (0 > D1 && D1 > D0)
        {
            ++flag;
        }
        if (0 > D2 && D2 > D0)
        {
            ++flag;
        }
    }
    if (flag == 2)
    {
        sa = D1 / D0;
        sb = D2 / D0;
    }

    // Out-line distance
    if (flag < 2)
    {
        // a11 >= 0, a22 >= 0
        s1 = 0.0;
        n2 = b2; // s2 = (b2 - a12 * s1 / a22
        flag = 0;
        if (n2 <= 0)
        {
            s2 = 0;
            ++flag;
        }
        if (n2 >= a22)
        {
            s2 = 1;
            ++flag;
        }
        if (flag == 0)
        {
            s2 = n2 / a22;
        }
        vecx(p1, g1, s1, x1);
        vecx(p2, g2, s2, x2);
        disa(x1, x2, da, & eq);
        sa = s1; sb = s2; dq = eq;
        s1 = 1.0;
        n2 = b2 - a12; // s2 = (b2 - a12 * s1) / a22
        flag = 0;
        if (n2 <= 0)
        {
            s2 = 0;
            ++flag;
        }
        if (n2 >= a22)
        {
            s2 = 1;
            ++flag;
        }
        if (flag == 0)
        {
            s2 = n2 / a22;
        }
        vecx(p1, g1, s1, x1);
        vecx(p2, g2, s2, x2);
        disa(x1, x2, da, & eq);

        if (eq < dq)
        {
            sa = s1; sb = s2; dq = eq;
        }
        s2 = 0;
        n1 = b1; // s1 = (b1 - a12 * s2) / a11
        flag = 0;
        if (n1 <= 0)
        {
            s1 = 0;
            ++flag;
        }
        if (n1 >= a11)
        {
            s1 = 1;
            ++flag;
        }
        if (flag == 0)
        {
            s1 = n1 / a11;
        }
        vecx(p1, g1, s1, x1);
        vecx(p2, g2, s2, x2);
        disa(x1, x2, da, &eq);
        if (eq < dq)
        {
            sa = s1; sb = s2; dq = eq;
        }
        s2 = 1.0;
        n1 = b1 - a12; // s1 = (b1 - a12 * s2) / a11
        flag = 0;
        if (n1 <= 0)
        {
            s1 = 0;
            ++flag;
        }
        if (n1 >= a11)
        {
            s1 = 1;
            ++flag;
        }
        if (flag == 0)
        {
            s1 = n1 / a11;
        }
        vecx(p1, g1, s1, x1);
        vecx(p2, g2, s2, x2);
        disa(x1, x2, da, & eq);
        if (eq < dq)
        {
            sa = s1; sb = s2; dq = eq;
        }
    }
    // For all sa = s1, sb = s2, refresh
    vecx(p1, g1, sa, x1);
    vecx(p2, g2, sb, x2);
    disa(x1, x2, da, & dq);

    gsl_vector_free(p1);
    gsl_vector_free(q1);
    gsl_vector_free(p2);
    gsl_vector_free(q2);

    gsl_vector_free(g1);
    gsl_vector_free(g2);
    gsl_vector_free(dp);
    gsl_vector_free(da);

    gsl_vector_free(x1);
    gsl_vector_free(x2);

    return sqrt(dq);
}


void * calc_distances(void * args)
{
    calc_input * input = static_cast<calc_input *>(args);
    input_stat * sample1 = input->sample1_beg;
    input_stat * sample2 = input->sample2_beg;
    output_stat * out_ptr = input->out_start;

    size_t num_bases = 4;
    size_t locus_id = input->locus_id;
    while (sample1 != input->sample1_end)
    {
            
        (*out_ptr).sample1_id = input->sample1_id;
        (*out_ptr).sample2_id = input->sample2_id;
        (*out_ptr).locus_id = locus_id;
        if ((*sample1).rank == 5 || (*sample2).rank == 5)
        {
            // one or both of the input combinations are unknown.
            // set metrics to dummy metrics
            (*out_ptr).dist_mode = -1.0;
            (*out_ptr).dist_lo = 0.0;
            (*out_ptr).dist_hi = 1000.0;
        }
        else
        {
            (*out_ptr).dist_mode =
                sqrt(
                     ssd(
                         (*sample1).mode, (*(sample1 + 1)).mode, (*(sample1 + 2)).mode, (*(sample1 + 3)).mode,
                         (*sample2).mode, (*(sample2 + 1)).mode, (*(sample2 + 2)).mode, (*(sample2 + 3)).mode)
                     );
                
            // construct two line segments in 3D.  Each segment
            // represents the constrained distribution of this locus'
            // composition for the sample.  The constraint assumes
            // that the two minor components c3 and c4 are fixed at
            // their mode, and that c1, the major component is
            // constrained to go from lo to min(hi, 1 - c3 - c4).  c2,
            // the second-highest component, is constrained by the
            // normalization, equals 1 - c1 - c3 - c4.  So, the line
            // segment thus constructed, expressed in 4D, and
            // expressed in component order, is:

            // {c1_lo, n, c3_m, c4_m} to
            // {min(c1_hi, 1 - c3_m - c4_m), n, c3_m, c4_m}


            float s1[8];
            float s2[8];

            float * p1, * p2;
            char const* bases = "ACGT";
            input_stat * stat;

            for (size_t s = 0; s != 2; ++s)
            {
                stat = (s == 0) ? sample1 : sample2;
                p1 = (s == 0) ? s1 : s2;
                p2 = (s == 0) ? s1 + 4 : s2 + 4;
                size_t r1 = 0, r2 = 0;

                for (size_t c = 0; c != 4; ++c)
                {
                    size_t b = strchr(bases, stat[c].base) - bases;
                    switch (stat[c].rank)
                    {
                    case 0: r1 = b; p1[b] = stat[c].lo; p2[b] = stat[c].hi; break;
                    case 1: r2 = b; p1[b] = 0; p2[b] = 0; break;
                    case 2: 
                    case 3: p1[b] = stat[c].mode; p2[b] = stat[c].mode; break;
                    }
                }
                p1[r2] = 1.0 - (p1[0] + p1[1] + p1[2] + p1[3]);
                p2[r1] -= (p2[0] + p2[1] + p2[2] + p2[3]) - 1.0; // adjust first component
                p2[r2] = 1.0 - (p2[0] + p2[1] + p2[2] + p2[3]);
            }
            
            // now calculate min distances between sample1 and sample2 segments
            (*out_ptr).dist_lo = nearest_dist(s1, s1 + 4, s2, s2 + 4);
            (*out_ptr).dist_hi = 10000.0;
        }   

        sample1 += num_bases;
        sample2 += num_bases;
        ++locus_id;
        ++out_ptr;
    }
    pthread_exit((void *) 0);
}


struct line_range
{
    size_t sample_input_size;
    std::vector<char *>::iterator start;
    std::vector<char *>::iterator end;
    input_stat * out;
};


void * init_input(void * args)
{
    line_range * input = static_cast<line_range *>(args);
    std::vector<char *>::iterator cur_line = input->start;
    size_t sample_input_size = input->sample_input_size;

    input_stat * stat;
    size_t sample_id;
    size_t locus_id;
    char hyp_base;
    size_t rank;
    float mode, lo, hi;
    size_t num_bases = 4;

    while (cur_line != input->end)
    {
        sscanf(*cur_line, "%zi\t%zi\t%c\t%zi\t%f\t%f\t%f",
               & locus_id,
               & sample_id,
               & hyp_base,
               & rank,
               & mode,
               & lo,
               & hi);

        stat = input->out + (sample_input_size * sample_id + (num_bases * locus_id));
        switch (hyp_base)
        {
        case 'A': break;
        case 'C': stat += 1; break;
        case 'G': stat += 2; break;
        case 'T': stat += 3; break;
        }
        (*stat).base = hyp_base;
        (*stat).rank = rank;
        (*stat).mode = mode;
        (*stat).lo = lo;
        (*stat).hi = hi;

        ++cur_line;
    }
    pthread_exit((void *) 0);
}



int main(int argc, char ** argv)
{
    size_t num_threads = static_cast<size_t>(atof(argv[1]));
    size_t num_samples = static_cast<size_t>(atof(argv[2]));
    size_t num_loci = static_cast<size_t>(atof(argv[3]));
    size_t chunk_size = static_cast<size_t>(atof(argv[4]));
    char * data_file = argv[5];
    char * comp_file = argv[6];

    size_t num_bases = 4;

    size_t total_input_size = num_samples * num_loci * num_bases;
    size_t sample_input_size = total_input_size / num_samples;

    input_stat * input_stats = new input_stat[total_input_size];
    input_stat * in_ptr = input_stats;
    input_stat * in_end = input_stats + total_input_size;

    size_t num_comparisons = 0;
    comparison * comparisons = new comparison[1000];
    comparison * comp_ptr = comparisons;
    // parse comparisons file
    FILE * comp_fh = fopen(comp_file, "r");
    if (comp_fh == NULL)
    {
        fprintf(stderr, "Error, couldn't open file %s\n", comp_file);
        exit(1);
    }

    while (! feof(comp_fh))
    {
        fscanf(comp_fh, "%zu\t%zu\n", & (*comp_ptr).sample1_id, & (*comp_ptr).sample2_id);
        comp_ptr++;
        num_comparisons++;
    }
    fclose(comp_fh);


    output_stat * output_stats = new output_stat[num_loci];
    output_stat * out_ptr = output_stats;
    output_stat * out_end = output_stats + num_loci;
    
    // since we are doing 'perfect insertion sort', we don't know which values will be missing.
    // therefore, initialize everything to the default of 'missing'
    while (in_ptr != in_end)
    {
        (*in_ptr++).rank = 5;
    }

    // load input
    FILE * in_fh = fopen(data_file, "r");
    if (in_fh == NULL)
    {
        fprintf(stderr, "Error, couldn't open data file %s\n", data_file);
        exit(1);
    }
    size_t max_line_length = 100;
    size_t nbytes_read;
    size_t nbytes_want = chunk_size - max_line_length;
    char * chunk_buffer = new char[chunk_size];
    char * dummy; // needed for 'find_complete_lines_nullify'

    while (! feof(in_fh))
    {
        nbytes_read = fread(chunk_buffer, 1, nbytes_want, in_fh);
        chunk_buffer[nbytes_read] = '\0';

        if (nbytes_read == nbytes_want)
        {
            // keep reading until we hit the next newline
            char * cur = chunk_buffer + nbytes_read;
            char * test = fgets(cur, max_line_length, in_fh);
            if (test != cur)
            {
                fprintf(stderr, "Error: unexpected line size\n");
                exit(2);
            }
            nbytes_read += strlen(cur);
        }
        else
        {
            // didn't get the requested bytes.  should be at the end of the file...
            if (chunk_buffer[nbytes_read - 1] != '\n')
            {
                fprintf(stderr, "Error: File doesn't end in a newline\n");
                exit(1);
            }
            else
            {
                // all good.  do nothing
            }
        }

        std::vector<char *> lines =
            FileUtils::find_complete_lines_nullify(chunk_buffer, & dummy);

        line_range * line_ranges = new line_range[num_threads];
        size_t work_load = lines.size() / num_threads;
        for (size_t t = 0; t != num_threads; ++t)
        {
            line_ranges[t].sample_input_size = sample_input_size;
            line_ranges[t].start = lines.begin() + (t * work_load);
            line_ranges[t].end = (t == (num_threads - 1)) ? lines.end() : lines.begin() + ((t + 1) * work_load);
            line_ranges[t].out = input_stats;
        }

        pthread_t * threads = new pthread_t[num_threads];
        for (size_t t = 0; t != num_threads; ++t)
        {
            int rc = pthread_create(&threads[t], NULL,
                                    & init_input,
                                    static_cast<void *>(& line_ranges[t]));
            assert(rc == 0);
        }

        for (size_t t = 0; t < num_threads; ++t) {
            int rc = pthread_join(threads[t], NULL);
            assert(0 == rc);
        }

        delete[] threads;
        delete[] line_ranges;

    }

    delete chunk_buffer;
    fclose(in_fh);

    // compute output stats
    input_stat * sample1;
    input_stat * sample2;
    input_stat * sample1_end;

    size_t max_out_line_size = 50;
    char * chunk_buffer_out = new char[max_out_line_size * num_loci];

    for (size_t c = 0; c != num_comparisons; ++c)
    {
        size_t sample1_id = comparisons[c].sample1_id;
        size_t sample2_id = comparisons[c].sample2_id;
        sample1 = input_stats + (sample_input_size * sample1_id);
        sample2 = input_stats + (sample_input_size * sample2_id);
        sample1_end = sample1 + sample_input_size;
        out_ptr = output_stats;

        // here, instead of a single loop, num_threads loops
        calc_input * calc_inputs = new calc_input[num_threads];
        size_t work_load = sample_input_size / (num_threads * num_bases);
        for (size_t t = 0; t != num_threads; ++t)
        {
            calc_inputs[t].sample1_id = sample1_id;
            calc_inputs[t].sample2_id = sample2_id;
            calc_inputs[t].locus_id = t * work_load;
            calc_inputs[t].sample1_beg = sample1 + (t * work_load * num_bases);
            calc_inputs[t].sample1_end = (t == (num_threads - 1)) ? sample1_end : sample1 + ((t + 1) * work_load * num_bases);
            calc_inputs[t].sample2_beg = sample2 + (t * work_load * num_bases);
            calc_inputs[t].out_start = out_ptr + (t * work_load);
        }

        pthread_t * threads = new pthread_t[num_threads];
        for (size_t t = 0; t != num_threads; ++t)
        {
            int rc = pthread_create(&threads[t], NULL,
                                    & calc_distances,
                                    static_cast<void *>(& calc_inputs[t]));
            assert(rc == 0);
        }

        for (size_t t = 0; t < num_threads; ++t) {
            int rc = pthread_join(threads[t], NULL);
            assert(rc == 0);
        }

        delete[] threads;
        delete[] calc_inputs;

        // write output stats
        out_ptr = output_stats;
        char * out_buf_ptr = chunk_buffer_out;
        while (out_ptr != out_end)
        {
            out_buf_ptr += 
            sprintf(out_buf_ptr,
                    "%u\t%u\t%u\t%6.5f\t%6.5f\t%6.5f\n",
                    (*out_ptr).sample1_id,
                    (*out_ptr).sample2_id,
                    (*out_ptr).locus_id,
                    (*out_ptr).dist_mode,
                    (*out_ptr).dist_lo,
                    (*out_ptr).dist_hi);
            ++out_ptr;
        }
        fwrite(chunk_buffer_out, 1, (out_buf_ptr - chunk_buffer_out), stdout);
        fflush(stdout);
        
    }
    
    delete chunk_buffer_out;
    delete input_stats;
    delete output_stats;
    delete comparisons;
}
