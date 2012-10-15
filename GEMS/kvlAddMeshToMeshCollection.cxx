/**
 * @file  kvlAddMeshToMeshCollection.cxx
 * @brief REPLACE_WITH_ONE_LINE_SHORT_DESCRIPTION
 *
 * REPLACE_WITH_LONG_DESCRIPTION_OR_REFERENCE
 */
/*
 * Original Author: Koen Van Leemput
 * CVS Revision Info:
 *    $Author: nicks $
 *    $Date: 2012/10/15 21:17:38 $
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


