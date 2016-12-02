#include "kvlAtlasMeshCollection.h"


int main( int argc, char** argv )
{
  // Sanity check on input
  if ( argc != 2 )
    {
    std::cerr << "Usage: "<< argv[ 0 ] << " inputMeshCollection" << std::endl;
    return -1;
    }


  // Read mesh collection from file
  kvl::AtlasMeshCollection::Pointer  inputCollection =  kvl::AtlasMeshCollection::New();
  inputCollection->Read( argv[ 1 ] );

  // Upsample
  kvl::AtlasMeshCollection::Pointer  outputCollection = inputCollection->GetUpsampled();

  // Write the result out
  outputCollection->Write( "upsampled.txt" );

  return 0;
};


