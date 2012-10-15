/**
 * @file  kvlViewMeshCollectionWithGUI.cxx
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
#include "kvlAtlasMeshViewingConsole.h"
#include "itkMGHImageIOFactory.h"


int main( int argc, char* argv[] )
{
  // Sanity check on input
  if ( argc < 2 )
  {
    std::cerr << argv[0] << " meshFileName [ sizeX sizeY sizeZ backgroundImage ]" << std::endl;
    return -1;
  }

  // Add support for MGH file format to ITK. An alternative way to add this by default would be
  // to edit ITK's itkImageIOFactory.cxx and explicitly adding it in the code there.
  itk::ObjectFactoryBase::RegisterFactory( itk::MGHImageIOFactory::New() );

  //
  kvl::AtlasMeshViewingConsole  console;
  if ( argc > 4 )
  {
    int size[ 3 ];
    for ( int i = 0; i < 3; i++ )
    {
      std::istringstream  sizeStream( argv[ i + 2 ] );
      sizeStream >> size[ i ];
    }

    std::string  backgroundImageFileName;
    if ( argc > 5 )
    {
      backgroundImageFileName = argv[ 5 ];
    }

    console.LoadMeshCollection( argv[ 1 ], size, backgroundImageFileName );
  }
  else
  {
    console.LoadMeshCollection( argv[ 1 ] );
  }

  console.Show();

  return 0;
};

