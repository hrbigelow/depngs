
// WARNING: this code overwrites the files reads.fastb, reads.qltout, and
// ViralTest.dot.

/*
  #ifndef FORCE_DEBUG
  #define NDEBUG
  #endif
*/

#include "Basevector.h"
#include "FetchReads.h"
#include "Feudal.h"
#include "MainTools.h"
#include "graph/Digraph.h"
#include "lookup/LookAlign.h"
#include "paths/HyperKmerPath.h"
#include "paths/KmerBaseBroker.h"
#include "paths/KmerPath.h"
#include "paths/ReadsToPathsCoreX.h"
#include "paths/Unipath.h"

int main(int argc, char ** argv )
{
  RunTime( );
  
  // Define the directory for the HIV data.
  
  String dir = argv[0];
  //"/seq/assembly_analysis/projects/HIV/20090324_XLR_Demi/"
  //"012251I/AssembleViral454_run";

  // Load the reads and quality scores.

  vecbasevector bases;
  vecqualvector quals;
  FetchReads( bases, 0, dir + "/fasta/reads.fasta" );
  ReadFastaQuals( dir + "/qual/reads.qual", quals );

  // Align the reads.  Here is where files are overwritten.

  bases.WriteAll( "reads.fastb" );
  cout << Date( ) << ": aligning" << endl;
  SystemSucceed( "QueryLookupTable K=12 MM=12 MC=0.15 SEQS=reads.fastb "
                 "L=/seq/assembly_analysis/projects/HIV/references/"
                 "HIV_K03455.fasta.lookuptable.lookup PARSEABLE=True > reads.qltout" );
  cout << Date( ) << ": done" << endl;
  vec<look_align> aligns;
  LoadLookAligns( "reads.qltout", aligns );

  // Filter out reads that are not on a small interval, so for the moment
  // we have less to look at.

  vec<Bool> to_remove( bases.size( ), True );
  for ( int i = 0; i < aligns.isize( ); i++ )
    {    const look_align& la = aligns[i];
      if ( la.pos2( ) >= 4000 && la.Pos2( ) <= 4500 )
        to_remove[la.query_id] = False;    }
  bases.EraseIf(to_remove);
  quals.EraseIf(to_remove);

  // ............................................................................

  int K = 20;
  int genome_size = 10000;

  bases.WriteAll( "xreads.fastb" );
  quals.WriteAll( "xreads.qualb" );
  vec<Bool> is_unique( bases.size( ), True );
  BinaryWrite3( "xreads.is_strong", is_unique );
  SystemSucceed( "FindErrors IN_HEAD=xreads OUT_HEAD=yreads WORKDIR=tmp DELETE=True" );
  bases.clear( );
  bases.ReadAll( "yreads.fastb" );
     
  vecbasevector ref( "/seq/assembly_analysis/projects/HIV/references/"
                     "HIV_K03455.fasta.lookuptable.fastb" );
  basevector reflet;
  reflet.SetToSubOf( ref[0], 3000, 2000 );
  bases.push_back(reflet);

  vecKmerPath paths, pathsrc, unipaths;
  vec<tagged_rpint> pathsdb, unipathsdb;
  ReadsToPathsCoreY( bases, K, paths, pathsrc, pathsdb );
  Unipath( paths, pathsrc, pathsdb, unipaths, unipathsdb );
  digraph A;
  BuildUnipathAdjacencyGraph( paths, pathsrc, pathsdb, unipaths, unipathsdb, A );
  HyperKmerPath h;
  BuildUnipathAdjacencyHyperKmerPath( K, A, unipaths, h );
  KmerBaseBroker kbb( K, paths, pathsrc, pathsdb, bases );

  // Shave off edges that go nowhere.

  // Build the graph.

  vec<Bool> looper( h.N( ), False );
  for ( int v = 0; v < h.N( ); v++ )
    looper[v] = h.LoopAt(v);

  // Go through two passes, for reverse and forward directions.

  vec<Bool> edges_to_trim( h.EdgeObjectCount( ), False );
  for ( int xpass = 1; xpass <= 2; xpass++ )
    {    h.Reverse( );
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
            goes_to_looper[v] = True;

          for ( int i = 0; i < paths.isize( ); i++ )
            {    int len = 0;
              for ( int j = 0; j < paths[i].isize( ); j++ )
                {    int x = paths[i][j];
                  if ( looper[ paths[i][j] ] ) goes_to_looper[v] = True;
                  if ( j == paths[i].isize( ) - 1 ) continue;
                  int y = paths[i][j+1];
                  int M = 0;
                  for ( int u = 0; u < h.From(x).isize( ); u++ )
                    {    if ( h.From(x)[u] != y ) continue;
                      const KmerPath& p = h.EdgeObjectByIndexFrom(x, u);
                      M = Max( M, p.KmerCount( ) );    }
                  len += M;    }
              extent_fw[v] = Max( extent_fw[v], len );    }    }

      // Identify edges to delete.
     
      vec<int> to_left, to_right;
      h.ToLeft(to_left), h.ToRight(to_right);
      for ( int v = 0; v < h.N( ); v++ )
        {    int biggest = 0;
          for ( int pass = 1; pass <= 2; pass++ )
            {    for ( int i = 0; i < h.From(v).isize( ); i++ )
                {    int w = h.From(v)[i];
                  int e = h.EdgeObjectIndexByIndexFrom( v, i );
                  int len = h.EdgeObjectByIndexFrom( v, i ).KmerCount( )
                    + extent_fw[w];
                  biggest = Max( biggest, len );
                  if ( pass == 1 ) continue;
                  if ( !goes_to_looper[w] 
                       && 10 * len < biggest && len < 100 )
                    {    cout << "kill edge " << e << endl;    
                      vec<int> todie;
                      todie.push_back(e);
                      while( todie.nonempty( ) )
                        {    int f = todie.back( );
                          todie.resize( todie.isize( ) - 1 );
                          if ( edges_to_trim[f] ) continue;
                          cout << "and kill edge " << f << endl;
                          edges_to_trim[f] = True;
                          int r = to_right[f];
                          for ( int z = 0; z < h.From(r).isize(); z++ )
                            {    todie.push_back(
                                                 h.EdgeObjectIndexByIndexFrom( 
                                                                              r, z ) );
                            }    }    }    }    }    }    }

  // Delete the edges.

  vec<int> to_delete;
  for ( int i = 0; i < edges_to_trim.isize( ); i++ )
    if ( edges_to_trim[i] ) to_delete.push_back(i);
  h.DeleteEdges(to_delete);
  h.RemoveDeadEdgeObjects( );

  // Rebuild the unipaths.

  vecbasevector all;
  vec<int> to_left, to_right;
  h.ToLeft(to_left), h.ToRight(to_right);
  for ( int e = 0; e < h.EdgeObjectCount( ); e++ )
    {    all.push_back_reserve( kbb.Seq( h.EdgeObject(e) ) );
      int v = to_right[e];
      for ( int i = 0; i < h.From(v).isize( ); i++ )
        {    basevector p1 = kbb.Seq( h.EdgeObject(e) );
          basevector p2 = kbb.Seq( h.EdgeObjectByIndexFrom( v, i ) );
          p1.resize( p1.isize( ) - (K-1) );
          all.push_back_reserve( Cat( p1, p2 ) );    }    }
  vecKmerPath newpaths, newpathsrc, newunipaths;
  vec<tagged_rpint> newpathsdb, newunipathsdb;
  ReadsToPathsCoreY( all, K, newpaths, newpathsrc, newpathsdb );
  Unipath( newpaths, newpathsrc, newpathsdb, newunipaths, newunipathsdb );
  KmerBaseBroker newkbb( K, newpaths, newpathsrc, newpathsdb, all );
  newkbb.Seq( newunipaths[113] ).Print( cout, 113 );
  newkbb.Seq( newunipaths[436] ).Print( cout, 436 );
  newkbb.Seq( newunipaths[146] ).Print( cout, 146 );
  newkbb.Seq( newunipaths[114] ).Print( cout, 114 );

  digraph newA;
  BuildUnipathAdjacencyGraph( newpaths, newpathsrc, newpathsdb, newunipaths, 
                              newunipathsdb, newA );
  HyperKmerPath newh;
  BuildUnipathAdjacencyHyperKmerPath( K, newA, newunipaths, newh );

  Ofstream( out, "ViralTest.dot" );
  newh.PrintSummaryDOT0w( out, True, False, True, 0, False );    }
