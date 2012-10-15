/**
 * @file  kvlUpsampleMesh.cxx
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


