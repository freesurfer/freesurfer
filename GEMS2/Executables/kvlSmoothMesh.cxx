#include "kvlAtlasMeshSmoother.h"


int main( int argc, char** argv )
{
  // Sanity check on input
  if ( argc != 3 )
    {
    std::cerr << "Usage: "<< argv[ 0 ] << " inputMeshCollection sigma" << std::endl;
    return -1;
    }


  // Read mesh collection from file
  kvl::AtlasMeshCollection::Pointer  inputCollection =  kvl::AtlasMeshCollection::New();
  if ( !inputCollection->Read( argv[ 1 ] ) )
    {
    std::cout << "Couldn't read mesh collection from file " << argv[ 1 ] << std::endl;
    exit( -1 );
    }

  // Read the input sigma
  std::istringstream  sigmaStream( argv[ 2 ] );
  float  sigma;
  sigmaStream >> sigma;
  
  // Smooth
  kvl::AtlasMeshSmoother::Pointer  smoother = kvl::AtlasMeshSmoother::New();
  smoother->SetMeshCollection( inputCollection );
  smoother->SetSigma( sigma );
  kvl::AtlasMeshCollection::Pointer  outputCollection = smoother->GetSmoothedMeshCollection();

  // Write the result out
  std::ostringstream  outputFileNameStream;
  outputFileNameStream << argv[ 1 ] << "_smoothed" << sigma << ".txt";
  if ( !outputCollection->Write( outputFileNameStream.str().c_str() ) )
    {
    std::cerr << "Couldn't write mesh collection to file " << outputFileNameStream.str() << std::endl;
    exit( -1 );
    }
    std::cout << "Just wrote mesh collection to file " << outputFileNameStream.str() << std::endl;
  

  return 0;
};


