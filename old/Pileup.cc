#include <set>
#include <string>
#include <sstream>
#include <vector>

#include "Basevector.h"
#include "pairwise_aligners/SmithWatFree.h"


struct ConsensusTable {
  ConsensusTable() : insert_count(0), insert_sizes(), 
                     insert_sequences(), delete_count(0){
    nucleotide_count[0] = 0;
    nucleotide_count[1] = 0;
    nucleotide_count[2] = 0;
    nucleotide_count[3] = 0;

    quality_sum[0] = 0;
    quality_sum[1] = 0;
    quality_sum[2] = 0;
    quality_sum[3] = 0;
    
  }

  int nucleotide_count[4]; //number of reads with each of these nucleotides aligned here
  int quality_sum[4];
  int insert_count; //number of reads with an insert just after this nucleotide
  std::multiset<int> insert_sizes;
  std::multiset<std::string> insert_sequences;
  
  int delete_count; //number of reads missing the following nucleotide
};




int main(int argc, char ** argv){

  if (argc == 1){
    printf("Pileup <reference_file.f> <reads_file.fb> <quals_file.qb> "
           "<max_sw_score> <output_file>\n");
    exit(0);
  }

  char * reference_file = argv[1];
  char * reads_file = argv[2];
  char * quals_file = argv[3];
  unsigned int max_sw_score = atoi(argv[4]);
  char * output_file = argv[5];
  char * score_histo_file = argv[6];

    //now do the alignment
  alignment pairwise_alignment_fw, pairwise_alignment_rc;
  int best_location;
  align pairwise_alignment;
  bool forward_aligned;

  vecbasevector reference_sequences;

  reference_sequences.ReadAll(String(reference_file));

  vecbasevector reads;
  reads.ReadAll(String(reads_file));

  vecqualvector quals;
  quals.ReadAll(String(quals_file));

  
  //arbitrarily choose the first of the reference sequences
  basevector const& reference_sequence = reference_sequences[0];

  std::vector<ConsensusTable> built_consensus(reference_sequence.size());

  FILE * score_histo_stream = fopen(score_histo_file, "w");

  for (int r = 0; r != reads.size(); ++r){

    if (r % 100 == 0) {
      printf("Processing read %i\n", r);
      fflush(stdout);
    }

    unsigned int fw_score = 
      SmithWatFree(reads[r], reference_sequence, best_location, pairwise_alignment_fw);

    basevector read_rc;
    read_rc.ReverseComplement(reads[r]);

    qualvector qual_rev(quals[r]);
    qual_rev.ReverseMe();

    unsigned int rc_score =
      SmithWatFree(read_rc, reference_sequence, best_location, pairwise_alignment_rc);
    
    forward_aligned = fw_score < rc_score;

    unsigned int & chosen_score = forward_aligned ? fw_score : rc_score;

    fprintf(score_histo_stream, "%i\n", chosen_score);

    if (chosen_score > max_sw_score) continue;

    alignment & chosen_alignment = 
      forward_aligned ? pairwise_alignment_fw : pairwise_alignment_rc;

    basevector & chosen_read =
      forward_aligned ? reads[r] : read_rc;

    qualvector & chosen_qual =
      forward_aligned ? quals[r] : qual_rev;

    //assume 'first sequence' (PackAlign.h docs) is the read, or 'query'
    align chosen_align(chosen_alignment);

    int target_pos = chosen_align.StartOnTarget();
    int query_pos = chosen_align.StartOnQuery();

    //The total position range on the target should span from
    //the initial target position to the total size of all blocks and
    //positive gaps (gaps on the first sequence, which mean nucleotides on the
    //target
    for (int b = 0; b != chosen_align.Nblocks(); ++b) {

      //process current gap
      int gaplen = abs(chosen_align.Gaps()(b));
      bool gap_on_target = chosen_align.Gaps()(b) < 0;
      if (gap_on_target) {
        built_consensus[target_pos].insert_count++;
        built_consensus[target_pos].insert_sizes.insert(gaplen);
        basevector read_insert; read_insert.SetToSubOf(chosen_read, 
                                                       query_pos, gaplen);
        built_consensus[target_pos].insert_sequences.insert
          (std::string(read_insert.ToString().c_str()));

      }

      else for (int n = 0; n != gaplen; ++n)
             built_consensus[target_pos++].delete_count++;

      //process current alignment block
      int block_length = chosen_align.Lengths()(b);
      for (int n = 0; n != block_length; ++n){
        built_consensus[target_pos].nucleotide_count[chosen_read[query_pos]]++;
        built_consensus[target_pos].quality_sum[chosen_read[query_pos]] +=
          chosen_qual[query_pos];

        ++query_pos;
        ++target_pos;
      }

    }


  }

  fclose(score_histo_stream);

  FILE * output_stream = fopen(output_file, "w");

  //print out the built consensus
  fprintf(output_stream, 
          "position\tconsensus_base\tA_count\tC_count\tG_count\tT_count\t"
          "Num_read_inserts_here\tNum_read_gaps\t[Insert_size_list]\n");

  char const* base_letters = "ACGT";
  char const* different = "DIFFERENT";
  char const* blank = "";
  

  for (size_t n = 0; n != built_consensus.size(); ++n){
    std::ostringstream insert_sizes_list;
    
    std::copy(built_consensus[n].insert_sizes.begin(),
              built_consensus[n].insert_sizes.end(), 
              std::ostream_iterator<int>(insert_sizes_list, ","));

    std::ostringstream insert_sequences_list;
    std::copy(built_consensus[n].insert_sequences.begin(),
              built_consensus[n].insert_sequences.end(), 
              std::ostream_iterator<std::string>(insert_sequences_list, ","));

    int max_quality_sum = built_consensus[n].quality_sum[0];
    int max_quality_index = 0;
    for (int i=1; i != 4; ++i){
      if (max_quality_sum < built_consensus[n].quality_sum[i]){
        max_quality_sum = built_consensus[n].quality_sum[i];
        max_quality_index = i;
      }
    }

    char const* if_different = max_quality_index == reference_sequence[n] ?
      blank : different;

    int bold_indices[] = { 0, 0, 0, 0 };
    int base_colors[] = { 31, 32, 33, 34 };
    bold_indices[max_quality_index] = 1;

    fprintf(output_stream, "%i\t\033[%im%c\033[0m\t", 
            static_cast<int>(n),
            base_colors[reference_sequence[n]],
            base_letters[reference_sequence[n]]);

    for (int b = 0; b != 4; ++b){
      fprintf(output_stream, "\033[%i;%im%5i(%i) ", 
              base_colors[b], bold_indices[b],
              built_consensus[n].nucleotide_count[b],
              built_consensus[n].quality_sum[b]);
    }
    fprintf(output_stream, "\033[0mD[%i] I[%i][%s]%s\n",
            built_consensus[n].delete_count,
            built_consensus[n].insert_count,
            insert_sequences_list.str().c_str(),
            if_different);

  }

  fclose(output_stream);
  
}
