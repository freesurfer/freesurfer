#include "kvlAtlasMeshCollection.h"


int main( int argc, char** argv )
{
  // Sanity check on input
  if ( argc != 3 )
    {
    std::cerr << "Usage: "<< argv[ 0 ] << " meshCollection pointParametersMeshCollection" << std::endl;
    return -1;
    }


  // Read mesh collection from file
  kvl::AtlasMeshCollection::Pointer  collection =  kvl::AtlasMeshCollection::New();
  collection->Read( argv[ 1 ] );

  // Read pointParameteresMeshCollection from file
  kvl::AtlasMeshCollection::Pointer  sourceCollection =  kvl::AtlasMeshCollection::New();
  sourceCollection->Read( argv[ 2 ] );

  // Copy point parameters
  kvl::AtlasMesh::PointDataContainer::Iterator  targetIt = collection->GetPointParameters()->Begin();
  kvl::AtlasMesh::PointDataContainer::ConstIterator  sourceIt = sourceCollection->GetPointParameters()->Begin();
  for ( ; targetIt != collection->GetPointParameters()->End(); ++targetIt, ++sourceIt )
    {
    targetIt.Value() = sourceIt.Value();
    }


  // Write the result out
  std::ostringstream  fileNameStream;
  fileNameStream << argv[ 1 ] << "_changedPointParameters.txt";
  if ( !collection->Write( fileNameStream.str().c_str() ) )
    {
    std::cerr << "Couldn't write mesh collection to file " << fileNameStream.str() << std::endl;
    exit( -1 );
    }
  std::cout << "Just wrote mesh collection to file " << fileNameStream.str() << std::endl;


  return 0;
};
