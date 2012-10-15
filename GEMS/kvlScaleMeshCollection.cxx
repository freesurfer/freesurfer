/**
 * @file  kvlScaleMeshCollection.cxx
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


int main( int argc, char* argv[] )
{

  // Check input arguments
  if ( argc < 3 )
  {
    std::cerr << "Usage: " << argv[ 0 ] << " meshCollectionFileName scaleFactor [ scaleFactor scaleFactor ]" << std::endl;
    exit( -1 );
  }

  // Retrieve input parameters
  std::ostringstream  argumentOStream;
  for ( int i = 1; i < argc; i++ )
  {
    argumentOStream << argv[ i ] << " ";
  }
  std::istringstream  argumentIStream( argumentOStream.str() );
  std::string  meshCollectionFileName;
  float  scaleFactors[ 3 ];
  argumentIStream >> meshCollectionFileName >> scaleFactors[ 0 ];
  scaleFactors[ 1 ] = scaleFactors[ 0 ];
  scaleFactors[ 2 ] = scaleFactors[ 0 ];
  if ( argc > 3 )
  {
    argumentIStream >> scaleFactors[ 1 ];
  }
  if ( argc > 4 )
  {
    argumentIStream >> scaleFactors[ 2 ];
  }

  // Read the mesh collection
  kvl::AtlasMeshCollection::Pointer  collection = kvl::AtlasMeshCollection::New();
  if ( !collection->Read( meshCollectionFileName.c_str() ) )
  {
    std::cerr << "Couldn't read mesh collection from file " << meshCollectionFileName << std::endl;
    exit( -1 );
  }

  // Scale the reference position and all positions
  std::cout << "Having " << collection->GetNumberOfMeshes() << " meshes" << std::endl;
  for ( int meshNumber = -1; meshNumber < static_cast< int >( collection->GetNumberOfMeshes() ); meshNumber++ )
  {
    // Get a pointer to the position
    kvl::AtlasMesh::PointsContainer::Pointer  position = 0;
    if ( meshNumber == -1 )
    {
      position = collection->GetReferencePosition();
      std::cout << "Doing reference position" << std::endl;
    }
    else
    {
      std::cout << "Doing position " << meshNumber << std::endl;
      position = collection->GetPositions()[ meshNumber ];
    }

    // Loop over all points and scale them
    for ( kvl::AtlasMesh::PointsContainer::Iterator  it = position->Begin();
          it != position->End(); ++it )
    {
      //std::cout << it.Value() << " -> ";
      for ( int i = 0; i < 3; i++ )
      {
        it.Value()[ i ] *= scaleFactors [ i ];
      }
      //std::cout << it.Value() << std::endl;
    }

  }


  // Write out
  std::ostringstream  outputFileNameStream;
  outputFileNameStream << meshCollectionFileName << "_scaledByFactors_"
                       << scaleFactors[ 0 ] << "x"
                       << scaleFactors[ 1 ] << "x"
                       << scaleFactors[ 2 ] << ".txt";
  if ( !collection->Write( outputFileNameStream.str().c_str() ) )
  {
    std::cerr << "Couldn't write output to file " << outputFileNameStream.str().c_str() << std::endl;
    exit( -1 );
  }
  std::cout << "Just wrote output to file " << outputFileNameStream.str().c_str() << std::endl;

  return 0;
};

