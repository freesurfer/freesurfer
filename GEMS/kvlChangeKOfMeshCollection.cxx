/**
 * @file  kvlChangeKOfMeshCollection.cxx
 * @brief REPLACE_WITH_ONE_LINE_SHORT_DESCRIPTION
 *
 * REPLACE_WITH_LONG_DESCRIPTION_OR_REFERENCE
 */
/*
 * Original Author: Koen Van Leemput
 * CVS Revision Info:
 *    $Author: nicks $
 *    $Date: 2012/10/15 21:17:39 $
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
