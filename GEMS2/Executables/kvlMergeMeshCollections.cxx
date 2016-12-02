#include "kvlAtlasMeshCollection.h"


int main( int argc, char** argv )
{
  // Sanity check on input
  if ( argc < 3 )
    {
    std::cerr << "Usage: "<< argv[ 0 ] << " collection1 collection2 [ collection3 ... ]" << std::endl;
    return -1;
    }

  // Read original mesh collection from file
  kvl::AtlasMeshCollection::Pointer  collection =  kvl::AtlasMeshCollection::New();
  if ( !collection->Read( argv[ 1 ] ) )
    {
    std::cerr << "Couldn't read mesh collection from file " << argv[ 1 ] << std::endl;
    exit( -1 );
    }

  // Loop over all other mesh collections
  for ( int i = 2; i < argc; i++ )
    {
    // Read mesh collection
    kvl::AtlasMeshCollection::Pointer  collectionToCopyFrom =  kvl::AtlasMeshCollection::New();
    if ( !collectionToCopyFrom->Read( argv[ i ] ) )
      {
      std::cerr << "Couldn't read mesh collection from file " << argv[ i ] << std::endl;
      exit( -1 );
      }


    for ( unsigned int meshNumber = 0; meshNumber < collectionToCopyFrom->GetNumberOfMeshes(); meshNumber++ )
      {
      // Append
      ( collection->GetPositions() ).push_back( collectionToCopyFrom->GetPositions()[ meshNumber ] );
      }
    }


  // Write the result out
  const std::string  outputFileName = "merged.txt";
  if ( !collection->Write( outputFileName.c_str() ) )
    {
    std::cerr << "Couldn't not write mesh collection to file " << outputFileName << std::endl;
    exit( -1 );
    }
  std::cout << "Just wrote mesh collection to file " << outputFileName << std::endl;


  return 0;
};


