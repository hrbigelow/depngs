// WARNING: this code overwrites the files reads.fastb, reads.qltout, and
// ViralTest.dot.

/*
  #ifndef FORCE_DEBUG
  #define NDEBUG
  #endif
*/

#include <sstream>

#include "Basevector.h"
#include "FetchReads.h"
#include "FastaFileset.h"
#include "Feudal.h"
#include "MainTools.h"
#include "graph/Digraph.h"
#include "graph/DigraphTemplate.h"
#include "lookup/LookAlign.h"
#include "paths/HyperKmerPath.h"
#include "paths/KmerBaseBroker.h"
#include "paths/KmerPath.h"
#include "paths/ReadsToPathsCoreX.h"
#include "paths/Unipath.h"
#include "pairwise_aligners/SmithWatFree.h"



/* Enumerates all paths (sequences of *edges*) defined by the sequence
   of <vertex_indices> that correspond to <graph>
   As a temporary structure, generates a vec<Edge>, representing the
   currently enumerated edge path.  When it gets to the last vertex,
   applies (process_edge_path) to the edge path, using <data> as any
   required auxiliary data.

*/
template <typename Graph> void EnumerateEdgePaths(Graph const& graph,
                                                  vec<int> const& vertex_indices,
                                                  vec<int> & edge_path,
                                                  void * auxiliary_data,
                                                  size_t current_edge_pos,
                                                  void (process_edge_path
                                                        )(Graph const& graph,
                                                          vec<int> const& edge_path,
                                                          void * data)){
  
    if (current_edge_pos == vertex_indices.size() - 1){
        process_edge_path(graph, edge_path, auxiliary_data);
    } else {
    
        vec<int> edges = graph.EdgesBetween(vertex_indices[current_edge_pos],
                                            vertex_indices[current_edge_pos + 1]);

        for (size_t e = 0; e != edges.size(); ++e){

            if (edge_path.size() <= current_edge_pos) 
                edge_path.resize(current_edge_pos+1);

            edge_path[current_edge_pos] = edges[e];

            EnumerateEdgePaths(graph, vertex_indices, edge_path, auxiliary_data,
                               ++current_edge_pos, process_edge_path);
        }
    }


}


//
void PrintCollapsedEdgePathsLengths(FILE * outfile,
                                    HyperKmerPath const& graph,
                                    vec<int> const& vertex_indices){
  
    vec<int> current_edges;
    for (size_t v = 1; v != vertex_indices.size(); ++v){
        current_edges = graph.EdgesBetween(vertex_indices[v-1], vertex_indices[v]);
    
        if (current_edges.size() > 1) fprintf(outfile, "{");
        fprintf(outfile, "%4i[%i]", current_edges[0], 
                graph.EdgeLength(current_edges[0]));
        for (size_t e = 1; e != current_edges.size(); ++e)
            fprintf(outfile, " %4i[%i]", current_edges[e], 
                    graph.EdgeLength(current_edges[e]));
    
        if (current_edges.size() > 1) fprintf(outfile, "}");
    }
    fprintf(outfile, "\n");

}


//Print the edge path such that edges sharing the same source and target
//vertices are printed in brackets
void PrintCollapsedEdgePaths(FILE * outfile,
                             HyperKmerPath const& graph,
                             vec<int> const& vertex_indices){
  
    vec<int> current_edges;
    for (size_t v = 1; v != vertex_indices.size(); ++v){
        current_edges = graph.EdgesBetween(vertex_indices[v-1], vertex_indices[v]);
    
        if (current_edges.size() > 1) fprintf(outfile, "{");
        fprintf(outfile, "%4i", current_edges[0]);
        for (size_t e = 1; e != current_edges.size(); ++e)
            fprintf(outfile, " %4i", current_edges[e]);
    
        if (current_edges.size() > 1) fprintf(outfile, "}");
    }
    fprintf(outfile, "\n");
  
}


void PrintEdgePath(HyperKmerPath const& graph,
                   vec<int> const& edge_path, 
                   void * data){
  
    std::pair<KmerBaseBroker *, int> * bases_and_K = 
        static_cast<std::pair<KmerBaseBroker *, int> * >(data);
  
    KmerBaseBroker & kbb = *(bases_and_K->first);
    int K = bases_and_K->second;
  
    //print out the edge indices, in collapsed form, with length
    //edge_index[length]

    fprintf(stdout, "%4i[%i]", edge_path[0], graph.EdgeLength(edge_path[0]));
    for (size_t e = 1; e != edge_path.size(); ++e)
        fprintf(stdout, " %4i[%i]", edge_path[e], graph.EdgeLength(edge_path[e]));
    fprintf(stdout, "\n");

    //   SuperBaseVector edge_sequence;
    //   for (size_t e = 0; e != edge_path.size(); ++e){
    //     edge_sequence = kbb.ToSequence(*edge_path[e]);
    //     //print out the edge_sequence
    //     for (int b = 0; b != edge_sequence.size(); ++b){
    //       basevector const& bases = edge_sequence.Seq(b);
    //       for (size_t b2 = 0; b2 != bases.size() - K; ++b2)
    //         fprintf(stdout, "%c", bases.at(b2));

    //       //what to do about the gaps?  they are expressed as pairs of min, max size
      
    //     }
    //     fprintf(stdout, " ");
    //   }
    //   fprintf(stdout, "\n");
}



struct path_anchors {
    std::vector<int> start_indices;
    std::vector<int> end_indices;
};


int main(int argc, char ** argv ){
    //RunTime( );
  
    //fasta input file
    String fastb_file = argv[1];
    int K = atoi(argv[2]);
    int printed_component_min_size = atoi(argv[3]);
    char * unipath_fasta_file = argv[4];
    char * digraph_dot_file = argv[5];
    char * digraph_path_file = argv[6];
    char * complete_path_file = argv[7];

    char digraph_dot_file2[1000];
    strcpy(digraph_dot_file2, digraph_dot_file);
    strcat(digraph_dot_file2, "2");

    //String qualb_file = argv[2];

    // Load the reads and quality scores.
    vecbasevector reads;
    reads.ReadAll( fastb_file );

    //   vecqualvector quals;
    //   quals.ReadAll( qualb_file );

    int genome_size = 10000;

    // including a reference sequence
  
    //   vecbasevector ref;
    //   FastFetchReads(ref, 0, ref_genome);
    //reads.push_back(ref[0]);

    vecKmerPath paths, pathsrc, unipaths;
    vec<tagged_rpint> pathsdb, unipathsdb;
    ReadsToPathsCoreY( reads, K, paths, pathsrc, pathsdb );
    Unipath( paths, pathsrc, pathsdb, unipaths, unipathsdb );
    digraph A;
    BuildUnipathAdjacencyGraph( paths, pathsrc, pathsdb, unipaths, unipathsdb, A );
    HyperKmerPath h;
    BuildUnipathAdjacencyHyperKmerPath( K, A, unipaths, h );
    KmerBaseBroker kbb( K, paths, pathsrc, pathsdb, reads );

    // Shave off edges that go nowhere.

    // Build the graph.

    vec<Bool> looper( h.N( ), False );
    for ( int v = 0; v < h.N( ); v++ )
    {
        looper[v] = h.LoopAt(v);
    }

    // Go through two passes, for reverse and forward directions.

    vec<Bool> edges_to_trim( h.EdgeObjectCount( ), False );
    for ( int xpass = 1; xpass <= 2; xpass++)
    {    
        h.Reverse( );
        vec<Bool> goes_to_looper( h.N( ), False );
        vec<int> extent_fw( h.N( ), 0 );
    
        // For each vertex, see how far one can go forward, and see if one
        // encounters a cycles.

        for ( int v = 0; v < h.N( ); v++ )
        {    
            // Find all paths.  We treat failure as "goes_to_looper",
            // even though this doesn't really make sense.
      
            vec< vec<int> > paths;
            const int maxpaths = 100;
            const int maxpushes = 100;
            if ( !h.AllPaths( v, -1, paths, maxpaths, False, maxpushes ) )
            {
                goes_to_looper[v] = True;
            }
      
            for ( int i = 0; i < paths.isize( ); i++ ) 
            {    
                int len = 0;
                for ( int j = 0; j < paths[i].isize( ); j++ )
                {    
                    int x = paths[i][j];
                    if ( looper[ paths[i][j] ] ) 
                    {
                        goes_to_looper[v] = True;
                    }
                    if ( j == paths[i].isize( ) - 1 ) 
                    {
                        continue;
                    }
                    int y = paths[i][j+1];
                    int M = 0;
                    for ( int u = 0; u < h.From(x).isize( ); u++ )
                    {
                        if ( h.From(x)[u] != y ) 
                        {
                            continue;
                        }
                        const KmerPath& p = h.EdgeObjectByIndexFrom(x, u);
                        M = Max( M, p.KmerCount( ) );    
                    }
                    len += M;    
                }
                extent_fw[v] = Max( extent_fw[v], len );    
            }    
        }
    
        // Identify edges to delete.
    
        vec<int> to_left, to_right;
        h.ToLeft(to_left), h.ToRight(to_right);
        for ( int v = 0; v < h.N( ); v++ ) 
        {
            int biggest = 0;
            for ( int pass = 1; pass <= 2; pass++ )
            {
                for ( int i = 0; i < h.From(v).isize( ); i++ )
                {
                    int w = h.From(v)[i];
                    int e = h.EdgeObjectIndexByIndexFrom( v, i );
                    int len = h.EdgeObjectByIndexFrom( v, i ).KmerCount( )
                        + extent_fw[w];
                    biggest = Max( biggest, len );
                    if ( pass == 1 ) 
                    {
                        continue;
                    }
                    if ( !goes_to_looper[w] 
                         && 10 * len < biggest && len < 100 )
                    {    
                        //cerr << "kill edge " << e << endl;    
                        vec<int> todie;
                        todie.push_back(e);
                        while( todie.nonempty( ) )
                        {    
                            int f = todie.back( );
                            todie.resize( todie.isize( ) - 1 );
                            if ( edges_to_trim[f] ) 
                            {
                                continue;
                            }
                            //cerr << "and kill edge " << f << endl;
                            edges_to_trim[f] = True;
                            int r = to_right[f];
                            for ( int z = 0; z < h.From(r).isize(); z++ )
                            {    
                                todie.push_back(h.EdgeObjectIndexByIndexFrom(r, z ));
                            } 
                        }
                    }    
                }    
            }    
        }    
    }

    // Delete the edges.
    printf("deleting edges.\n");

    vec<int> to_delete;
    for ( int i = 0; i < edges_to_trim.isize( ); i++ )
        if ( edges_to_trim[i] ) to_delete.push_back(i);
    h.DeleteEdges(to_delete);
    h.RemoveDeadEdgeObjects( );

    // Rebuild the unipaths.
    printf("rebuilding unipaths.\n");

    vecbasevector all;
    vec<int> to_left, to_right;
    h.ToLeft(to_left), h.ToRight(to_right);
    for ( int e = 0; e < h.EdgeObjectCount( ); e++ ) 
    {
        all.push_back_reserve( kbb.Seq( h.EdgeObject(e) ) );
        int v = to_right[e];
        for ( int i = 0; i < h.From(v).isize( ); i++ )
        {
            basevector p1 = kbb.Seq( h.EdgeObject(e) );
            basevector p2 = kbb.Seq( h.EdgeObjectByIndexFrom( v, i ) );
            p1.resize( p1.isize( ) - (K-1) );
            all.push_back_reserve( Cat( p1, p2 ) );    
        }    
    }
    
    vecKmerPath newpaths, newpathsrc, newunipaths;
    vec<tagged_rpint> newpathsdb, newunipathsdb;
    ReadsToPathsCoreY( all, K, newpaths, newpathsrc, newpathsdb );
    Unipath( newpaths, newpathsrc, newpathsdb, newunipaths, newunipathsdb );
    KmerBaseBroker newkbb( K, newpaths, newpathsrc, newpathsdb, all );

    digraph newA;
    BuildUnipathAdjacencyGraph( newpaths, newpathsrc, newpathsdb, newunipaths, 
                                newunipathsdb, newA );
    HyperKmerPath unipath_in_kmers;
    BuildUnipathAdjacencyHyperKmerPath( K, newA, newunipaths, unipath_in_kmers );

    //find the start and end of the graph.
    ofstream digraph_outstream(digraph_dot_file);
    Bool label_contigs = True;
    Bool label_vertices = True;
    Bool label_edges = False;


    vec<double> edge_lengths(unipath_in_kmers.EdgeObjectCount(), 1.0);
    vec<vec<int> > edges;

    //these will just be the actual string representations of the
    vec<vec<String> > edge_labels;
    unipath_in_kmers.ComponentEdges(edges);

    edge_labels.resize(edges.size());
    for (size_t i = 0; i < edges.size(); ++i){
        edge_labels[i].resize(edges[i].size());
        for (size_t j = 0; j < edges[i].size(); ++j){
            edge_labels[i][j] = ToString(edges[i][j]);
        }
    }

    //Print the DOT file and fasta sequences from the hyper_kmerpath
    equiv_rel kmer_components;
    unipath_in_kmers.ComponentRelation(kmer_components);
    vec<int> reps;
    kmer_components.OrbitRepsAlt(reps);

    vec<int> ComponentsToPrint;
    for (size_t i = 0; i != reps.size(); ++i)
    {
        if (kmer_components.OrbitSize(reps[i]) >= printed_component_min_size)
        {
            ComponentsToPrint.push_back(reps[i]);
        }
    }
  
    unipath_in_kmers.PrintSummaryDOT0w( digraph_outstream, True, False, True, &ComponentsToPrint, False );    
  
    //unipath_in_kmers.DOT(digraph_outstream);
    digraph_outstream.close();


    //also output unipath_in_kmers in a simple edge-labeled DOT format.
    FILE * digraph_outstream2 = fopen(digraph_dot_file2, "w");
    fprintf(digraph_outstream2, "digraph G {\n\n");
    fprintf(digraph_outstream2, "rankdir=LR;\n");
    fprintf(digraph_outstream2, "node [width=0.1,height=0.1,fontsize=10,shape=point];\n");
    fprintf(digraph_outstream2, "edge [fontsize=12];\n");
    for (int v = 0; v < unipath_in_kmers.N(); ++v){
        vec<int> const& v_targets = unipath_in_kmers.From()[v];
        for ( int j = 0; j < v_targets.isize( ); j++ ){
            fprintf(digraph_outstream2, "%i -> %i [label=%i]\n", v, v_targets[j],
                    unipath_in_kmers.FromEdgeObj(v)[j]);
        }
    }
    fprintf(digraph_outstream2, "\n}\n");
    fclose(digraph_outstream2);


    std::string dot_topo_command =
        std::string("/bin/bash -c \"") + " dot -Tplain " + digraph_dot_file + 
        " -o >(grep node /dev/stdin | cut -f 2,3 -d ' ' > " +
        digraph_dot_file + ".topsort)\"";

    system(dot_topo_command.c_str());

    FILE * digraph_topo_sort = 
        fopen(std::string(std::string(digraph_dot_file) + ".topsort").c_str(), "r");

    int node_number;
    float node_x_position;
    std::multimap<float, int> partial_node_order;
    std::multimap<float, int>::iterator node_order_iter;
    typedef std::multimap<float, int>::iterator MITER;

    while (! feof(digraph_topo_sort)){
        fscanf(digraph_topo_sort, "%i %f\n", &node_number, &node_x_position);
        partial_node_order.insert(std::make_pair(node_x_position, node_number));
    }
    fclose(digraph_topo_sort);


    HyperBasevector new_hyper_basevector(unipath_in_kmers, newkbb);

    printf("Printing unipath fasta file %s.\n", unipath_fasta_file);
    FILE * unipath_fasta_stream = fopen(unipath_fasta_file, "w");
    for (int e = 0; e < new_hyper_basevector.EdgeObjectCount( ); e++)
        fprintf(unipath_fasta_stream, ">%i\n%s\n", e, 
                ToString(new_hyper_basevector.EdgeObject(e)).c_str() + K - 1);
  
    fclose(unipath_fasta_stream);

    //new_hyper_basevector.DumpFasta(unipath_fasta_file);


    //now print full sequences

    //   equiv_rel connected_components;
    //   newA.ComponentRelation(connected_components);

    float start_node_xpos = (*partial_node_order.begin()).first;
    float end_node_xpos = (*partial_node_order.rbegin()).first;
    printf("Got start position %f and end position %f.\n",
           start_node_xpos, end_node_xpos);

    std::pair<MITER, MITER> start_node_range, end_node_range;
  
    start_node_range = partial_node_order.equal_range(start_node_xpos);
    end_node_range = partial_node_order.equal_range(end_node_xpos);

    path_anchors anchors;

    int num_components = kmer_components.OrbitCount();
    std::map<int, path_anchors> component_anchors;
    std::map<int, path_anchors>::iterator anchors_iter;

    int component_number;
    int node_index;
    for (MITER start_iter = start_node_range.first;
         start_iter != start_node_range.second; ++start_iter){
        node_index = (*start_iter).second;
        component_number = kmer_components.ClassId(node_index);
        component_anchors[component_number].start_indices.push_back(node_index);
    }

    for (MITER end_iter = end_node_range.first;
         end_iter != end_node_range.second; ++end_iter){
        node_index = (*end_iter).second;
        component_number = kmer_components.ClassId(node_index);
        component_anchors[component_number].end_indices.push_back(node_index);
    }

    bool allow_self_loop = true;
    vec<vec<int> > complete_paths;
    int c;

    printf("Tracing all paths from each pair of start and end anchors\n");
    for (anchors_iter = component_anchors.begin();
         anchors_iter != component_anchors.end(); ++anchors_iter){
        path_anchors & anchors = (*anchors_iter).second;
        for (size_t sn = 0; sn != anchors.start_indices.size(); ++sn)
            for (size_t en = 0; en != anchors.end_indices.size(); ++en){
                vec<vec<int> > complete_path;
                int start_index = anchors.start_indices[sn];
                int end_index = anchors.end_indices[en];

                printf("Tracing all paths from start index %i to end index %i\n",
                       start_index, end_index);

                unipath_in_kmers.AllPaths(start_index, end_index, 
                                          complete_path, -1, (! allow_self_loop), -1);

                complete_paths.insert(complete_paths.end(), 
                                      complete_path.begin(), complete_path.end());
            }
        
    }    


    //build a basevector of the first complete path
    vec<int> current_edge;
    basevector current_side;

    printf("Printing %i complete paths to %s.\n", 
           static_cast<int>(complete_paths.size()),
           complete_path_file);

    FILE * complete_path_stream = fopen(complete_path_file, "w");
  
    for (size_t p = 0; p != complete_paths.size(); ++p){
        vec<int> & vertex_index = complete_paths[p];
        basevector complete_path;
        for (size_t v = 1; v != vertex_index.size(); ++v){
            current_edge = unipath_in_kmers.EdgesBetween(vertex_index[v-1], vertex_index[v]);
            current_side = new_hyper_basevector.EdgeObject(current_edge[0]);
            current_side.SetToSubOf(current_side, 0, current_side.size(),
                                    current_side.size() + ((v == 1) ? 0 : K-1));
            complete_path = Cat(complete_path, current_side);
      
        }
        fprintf(complete_path_stream, ">%i\n%s\n", static_cast<int>(p), 
                complete_path.ToString().c_str());
    
    }
    fclose(complete_path_stream);


    //print out the path components in fasta format
    FILE * digraph_path_stream = fopen(digraph_path_file, "w");
    for (size_t p = 0; p != complete_paths.size(); ++p){
        vec<int> & this_vertex_path = complete_paths[p];
        PrintCollapsedEdgePaths(digraph_path_stream, unipath_in_kmers, this_vertex_path);
    }
    fclose(digraph_path_stream);

    //the paths are given in vertex indices.  enumerate all combinations of
    //edges between them
    vec<int> current_edge_path;
  
    //std::pair<KmerBaseBroker *, int> bases_and_K = std::make_pair(&newkbb, K);

    //   for (size_t p = 0; p != complete_paths.size(); ++p){
    //     vec<int> & this_vertex_path = complete_paths[p];
    //     void * new_kbb_ptr = static_cast<void*>(&bases_and_K);
    
    //     EnumerateEdgePaths(unipath_in_kmers,
    //                                  this_vertex_path, current_edge_path,
    //                                  new_kbb_ptr, 
    //                                  0, &PrintEdgePath);
    //   }
    


}  
  
