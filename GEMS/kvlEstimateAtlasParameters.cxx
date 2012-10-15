/**
 * @file  kvlEstimateAtlasParameters.cxx
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
#include "kvlAtlasParameterEstimationConsole.h"



int main( int argc, char** argv )
{
  //
  if ( argc < 4 )
  {
    std::cerr << "Usage: " << argv[ 0 ] << " useGaussians ignoreLastImage fileName1 [ fileName2 ... ]" << std::endl;
    exit( -1 );
  }

  // Collect useGaussians
  std::istringstream  useGaussiansStream( argv[ 1 ] );
  bool  useGaussians;
  useGaussiansStream >> useGaussians;

  // Collect ignoreLastImage
  std::istringstream  ignoreLastImageStream( argv[ 2 ] );
  bool  ignoreLastImage;
  ignoreLastImageStream >> ignoreLastImage;


  // Collect input file names
  std::vector< std::string >  fileNames;
  for ( int i=3; i<argc; i++ )
  {
    fileNames.push_back( std::string( argv[ i ] ) );
  }

  //
  kvl::AtlasParameterEstimationConsole  console;
  console.SetLabelImages( fileNames, useGaussians, ignoreLastImage );
  console.Show();

  return 0;
};

