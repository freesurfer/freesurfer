/**
 * @file  kvlSegmentWithoutGUI.cxx
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
#include "kvlAtlasMeshSegmentationDriver.h"


int main( int argc, char** argv )
{

  // Sanity check on input
  if ( argc != 2 )
  {
    std::cerr << "Usage: " << argv[ 0 ] << " setUpFileName" << std::endl;
    return( -1 );
  }


  //
  try
  {

    kvl::AtlasMeshSegmentationDriver::Pointer  driver = kvl::AtlasMeshSegmentationDriver::New();
    driver->SetUp( argv[ 1 ] );

    driver->Segment();
  }
  catch( itk::ExceptionObject& e )
  {
    std::cerr << e << std::endl;
    return -1;
  }

  return 0;
};


