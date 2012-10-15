/**
 * @file  kvlValidateMeshCollection.cxx
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
#include "kvlAtlasMeshCollectionValidator.h"


int main( int argc, char* argv[] )
{

  if ( argc != 2 )
  {
    std::cerr << "Usage: " << argv[ 0 ] << " meshFileName" << std::endl;
    exit( -1 );
  }

  // Read the mesh collection
  kvl::AtlasMeshCollection::Pointer  collection = kvl::AtlasMeshCollection::New();
  collection->Read( argv[ 1 ] );

  // Validate collapsed edge mesh collection
  kvl::AtlasMeshCollectionValidator::Pointer  validator = kvl::AtlasMeshCollectionValidator::New();
  validator->Validate( collection );

  return 0;
};

