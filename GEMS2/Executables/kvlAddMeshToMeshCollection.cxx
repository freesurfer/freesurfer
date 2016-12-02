#include "kvlAtlasMeshCollection.h"


int main( int argc, char** argv )
{
  // Sanity check on input
  if ( argc < 3 )
    {
    std::cerr << "Usage: "<< argv[ 0 ] << " originalCollection collectionToCopyFrom [ meshNumber=0 ]" << std::endl;
    return -1;
    }


  // Read original mesh collection from file
  kvl::AtlasMeshCollection::Pointer  originalCollection =  kvl::AtlasMeshCollection::New();
  originalCollection->Read( argv[ 1 ] );

  // Read mesh collection to copy from from file
  kvl::AtlasMeshCollection::Pointer  collectionToCopyFrom =  kvl::AtlasMeshCollection::New();
  collectionToCopyFrom->Read( argv[ 2 ] );

  // Retrieve the mesh number input
  int meshNumber = 0;
  if ( argc > 3 )
    {
    std::istringstream  meshNumberStream( argv[ 3 ] );
    meshNumberStream >> meshNumber;
    }
  std::cout << "Doing meshNumber: " << meshNumber << std::endl;


  // Append
  ( originalCollection->GetPositions() ).push_back( collectionToCopyFrom->GetPositions()[ meshNumber ] );


  // Write the result out
  const std::string  outputFileName = "appended.txt";
  if ( !originalCollection->Write( outputFileName.c_str() ) )
    {
    std::cerr << "Couldn't not write mesh collection to file " << outputFileName << std::endl;
    exit( -1 );
    }
  std::cout << "Just wrote mesh collection to file " << outputFileName << std::endl;


  return 0;
};


