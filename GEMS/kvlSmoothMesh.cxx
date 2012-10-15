/**
 * @file  kvlSmoothMesh.cxx
 * @brief REPLACE_WITH_ONE_LINE_SHORT_DESCRIPTION
 *
 * REPLACE_WITH_LONG_DESCRIPTION_OR_REFERENCE
 */
/*
 * Original Author: Koen Van Leemput
 * CVS Revision Info:
 *    $Author: nicks $
 *    $Date: 2012/10/15 21:17:40 $
 *    $Revision: 1.3 $
 *
 * Copyright © 2011 The General Hospital Corporation (Boston, MA) "MGH"
 *
 * Terms and conditions for use, reproduction, distribution and contribution
 * are found in the 'FreeSurfer Software License Agreement' contained
 * in the file 'LICENSE' found in the FreeSurfer distribution, and here:
 *
 * https://surfer.nmr.mgh.harvard.edu/fswiki/FreeSurferSoftwareLicense
 *
 * Reporting: freesurfer@nmr.mgh.harvard.edu
 *
 */
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


