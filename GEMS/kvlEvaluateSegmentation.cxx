/**
 * @file  kvlEvaluateSegmentation.cxx
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
#include "kvlSegmentationEvaluationConsole.h"


int main( int argc, char* argv[] )
{
  // Sanity check on input
  if ( argc < 4 )
  {
    std::cerr << argv[0] << " imageFileName overlayImageFileName1 overlayImageFileName2 [ overlayLookupTable ]" << std::endl;
    return -1;
  }

  //
  kvl::SegmentationEvaluationConsole  console;
  console.LoadImage( argv[ 1 ] );
  console.LoadOverlayImage1( argv[ 2 ] );
  console.LoadOverlayImage2( argv[ 3 ] );

  if ( argc > 4 )
  {
    console.LoadOverlayLookupTable( argv[ 4 ] );
  }

  console.Show();

  return 0;
};

