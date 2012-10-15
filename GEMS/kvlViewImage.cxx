/**
 * @file  kvlViewImage.cxx
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
#include "kvlImageViewingConsole.h"
#include "itkMGHImageIOFactory.h"


int main( int argc, char* argv[] )
{
  // Sanity check on input
  if ( argc < 2 )
  {
    std::cerr << argv[0] << " imageFileName [ overlayImageFileName overlayLookupTable meshCollectionFileName meshNumber ]" << std::endl;
    return -1;
  }

  // Add support for MGH file format to ITK. An alternative way to add this by default would be
  // to edit ITK's itkImageIOFactory.cxx and explicitly adding it in the code there.
  itk::ObjectFactoryBase::RegisterFactory( itk::MGHImageIOFactory::New() );


  //
  kvl::ImageViewingConsole  console;
  console.LoadImage( argv[ 1 ] );

  if ( argc > 2 )
  {
    console.LoadOverlayImage( argv[ 2 ] );
  }

  if ( argc > 3 )
  {
    console.LoadOverlayLookupTable( argv[ 3 ] );
  }

  if ( argc > 4 )
  {
    int  meshNumber = -1;
    if ( argc > 5 )
    {
      std::istringstream  meshNumberStream( argv[ 5 ] );
      meshNumberStream >> meshNumber;
    }

    console.LoadMesh( argv[ 4 ], meshNumber );
  }

  console.Show();

  return 0;
};

