#include "kvlAtlasMeshCollection.h"


int main( int argc, char** argv )
{
  // Sanity check on input
  if ( argc != 3 )
    {
    std::cerr << "Usage: "<< argv[ 0 ] << " inputMeshCollection K" << std::endl;
    return -1;
    }


  // Read mesh collection from file
  kvl::AtlasMeshCollection::Pointer  collection =  kvl::AtlasMeshCollection::New();
  collection->Read( argv[ 1 ] );

  // Read the input K
  std::istringstream  KStream( argv[ 2 ] );
  float  K;
  KStream >> K;
  
  // Set K
  collection->SetK( K );

  // Write the result out
  std::ostringstream  fileNameStream;
  fileNameStream << argv[ 1 ] << "_changedKTo" << K << ".txt";
  if ( !collection->Write( fileNameStream.str().c_str() ) )
    {
    std::cerr << "Could not write mesh collection to " << fileNameStream.str() << std::endl;
    exit( -1 );
    }
  std::cout << "Just wrote mesh collection to " << fileNameStream.str() << std::endl;


  return 0;
};
