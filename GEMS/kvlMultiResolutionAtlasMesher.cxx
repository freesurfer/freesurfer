/**
 * @file  kvlMultiResolutionAtlasMesher.cxx
 * @brief REPLACE_WITH_ONE_LINE_SHORT_DESCRIPTION
 *
 * REPLACE_WITH_LONG_DESCRIPTION_OR_REFERENCE
 */
/*
 * Original Author: Koen Van Leemput
 * CVS Revision Info:
 *    $Author: nicks $
 *    $Date: 2012/10/15 21:17:40 $
 *    $Revision: 1.4 $
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
#include "kvlMultiResolutionAtlasMesher.h"

#include <fstream>
#include "tetgen.h"


namespace kvl
{


//
//
//
MultiResolutionAtlasMesher
::MultiResolutionAtlasMesher()
{
  m_InitialSize[ 0 ] = 10;
  m_InitialSize[ 1 ] = 10;
  m_InitialSize[ 2 ] = 10;

  m_NumberOfUpsamplingSteps = 1;

  m_TryToBeSparse = false;

  m_InitialStiffnesses[ 0 ] = 0.1f;
  m_InitialStiffnesses[ 1 ] = 0.1f;
  m_InitialStiffnesses[ 2 ] = 0.1f;
  m_InitialStiffnesses[ 3 ] = 0.1f;
  m_InitialStiffnesses[ 4 ] = 0.1f;

  m_Estimator = AtlasParameterEstimator::New();

  m_NumberOfClasses = 0;
  m_NumberOfMeshes = 0;

  m_DomainSize[ 0 ] = 0;
  m_DomainSize[ 1 ] = 0;
  m_DomainSize[ 2 ] = 0;

  m_Current = 0;
  m_Hexahedra = 0;

}



//
//
//
MultiResolutionAtlasMesher
::~MultiResolutionAtlasMesher()
{
  // Clean up the hexahedral cells we've created
  for ( AtlasMesh::CellsContainer::Iterator  hexIt  = m_Hexahedra->Begin();
        hexIt != m_Hexahedra->End(); ++hexIt )
  {
    delete hexIt.Value();
  }

}



//
//
//
void
MultiResolutionAtlasMesher
::PrintSelf( std::ostream& os, itk::Indent indent ) const
{



}



//
//
//
void
MultiResolutionAtlasMesher
::SetLabelImages( const std::vector< LabelImageType::ConstPointer >& labelImages )
{
  m_LabelImages = labelImages;

}








//
//
//
void
MultiResolutionAtlasMesher
::SetUp( const unsigned int* size, const float* initialStiffnesses )
{
  m_InitialSize[ 0 ] = size[ 0 ];
  m_InitialSize[ 1 ] = size[ 1 ];
  m_InitialSize[ 2 ] = size[ 2 ];

  for ( unsigned int i = 0; i<5; i++ )
  {
    m_InitialStiffnesses[ i ] = initialStiffnesses[ i ];
  }

  this->SetUp();

}


//
//
//
void
MultiResolutionAtlasMesher
::SetUp()
{
  // Pass the label images onto the estimator
  m_Estimator->SetLabelImages( m_LabelImages );

  // Don't build a mesh if no label images have been set
  if ( m_Estimator->GetNumberOfLabelImages() == 0 )
  {
    return;
  }


  // Retrieve initial mesh parameters
  AtlasParameterEstimator::LabelImageType::SizeType  labelImageSize =
    m_Estimator->GetLabelImage( 0 )->GetLargestPossibleRegion().GetSize();
  m_DomainSize[ 0 ] = static_cast< unsigned int >( labelImageSize[ 0 ] );
  m_DomainSize[ 1 ] = static_cast< unsigned int >( labelImageSize[ 1 ] );
  m_DomainSize[ 2 ] = static_cast< unsigned int >( labelImageSize[ 2 ] );

  m_NumberOfClasses = m_Estimator->GetNumberOfClasses();

  m_NumberOfMeshes = m_Estimator->GetNumberOfLabelImages();


  // Set up hexahedra and reference position
  m_Hexahedra = AtlasMesh::CellsContainer::New();
  MeshSourceType::Pointer  meshSource = MeshSourceType::New();
  for ( unsigned int  x = 0; x < m_InitialSize[ 0 ]-1; x++ )
  {
    for ( unsigned int  y = 0; y < m_InitialSize[ 1 ]-1; y++ )
    {
      for ( unsigned int  z = 0; z < m_InitialSize[ 2 ]-1; z++ )
      {
        // Construct the eight corner coordinates
        float  x1 = static_cast< float >( x ) * static_cast< float >( m_DomainSize[ 0 ] - 1 )
                    / static_cast< float >( m_InitialSize[ 0 ] - 1 );
        float  y1 = static_cast< float >( y ) * static_cast< float >( m_DomainSize[ 1 ] - 1 )
                    / static_cast< float >( m_InitialSize[ 1 ] - 1 );
        float  z1 = static_cast< float >( z ) * static_cast< float >( m_DomainSize[ 2 ] - 1 )
                    / static_cast< float >( m_InitialSize[ 2 ] - 1 );

        float  x2 = static_cast< float >( x + 1 ) * static_cast< float >( m_DomainSize[ 0 ] - 1 )
                    / static_cast< float >( m_InitialSize[ 0 ] - 1 );
        float  y2 = static_cast< float >( y + 1 ) * static_cast< float >( m_DomainSize[ 1 ] - 1 )
                    / static_cast< float >( m_InitialSize[ 1 ] - 1 );
        float  z2 = static_cast< float >( z + 1 ) * static_cast< float >( m_DomainSize[ 2 ] - 1 )
                    / static_cast< float >( m_InitialSize[ 2 ] - 1 );

        const float  p0[] = { x1, y1, z1 };
        const float  p1[] = { x2, y1, z1 };
        const float  p2[] = { x1, y2, z1 };
        const float  p3[] = { x2, y2, z1 };
        const float  p4[] = { x1, y1, z2 };
        const float  p5[] = { x2, y1, z2 };
        const float  p6[] = { x1, y2, z2 };
        const float  p7[] = { x2, y2, z2 };

        this->AddHexahedron( meshSource, m_Hexahedra, p0, p1, p2, p3, p4, p5, p6, p7 );
      }

    }

  }
  AtlasMesh::PointsContainer::Pointer  referencePosition = meshSource->GetOutput()->GetPoints();


  // Make copies of the reference position
  std::vector< AtlasMesh::PointsContainer::Pointer >  positions;
  for ( unsigned int i = 0; i < m_NumberOfMeshes; i++ )
  {
    AtlasMesh::PointsContainer::Pointer  target = AtlasMesh::PointsContainer::New();

    for ( AtlasMesh::PointsContainer::ConstIterator  sourceIt = referencePosition->Begin();
          sourceIt != referencePosition->End(); ++sourceIt )
    {
      target->InsertElement( sourceIt.Index(), sourceIt.Value() );
    }

    // Push back
    positions.push_back( target );
  }


  // Create a mesh collection according to the reference position and the positions
  m_Current = this->GetMeshCollection( referencePosition, positions, m_InitialStiffnesses[ 0 ] );

  // Now initialize the estimator with the mesh
  m_Estimator->SetInitialMeshCollection( m_Current );


}




//
//
//
const MultiResolutionAtlasMesher::LabelImageType*
MultiResolutionAtlasMesher
::GetLabelImage( unsigned int labelImageNumber ) const
{
  // Sanity check
  if ( labelImageNumber >= m_LabelImages.size() )
  {
    return 0;
  }

  return m_LabelImages[ labelImageNumber ];
}




//
//
//
void
MultiResolutionAtlasMesher
::Go()
{

  //
  for ( unsigned int upsamplingStepNumber = 0; upsamplingStepNumber <= m_NumberOfUpsamplingSteps; upsamplingStepNumber++ )
  {
    std::cout << "running for upsamplingStepNumber: " << upsamplingStepNumber << std::endl;

    // Estimate
    std::cout << "       estimating..." << std::endl;
    m_Estimator->Estimate();

    // If this is not the final resolution yet, upsample mesh collection
    if ( upsamplingStepNumber != m_NumberOfUpsamplingSteps )
    {
      std::cout << "       upsampling..." << std::endl;
      this->Upsample();
      m_Current->SetK( m_InitialStiffnesses[ upsamplingStepNumber + 1 ] );
      m_Estimator->SetInitialMeshCollection( m_Current );
    }

  } // End loop over upsampling steps


}




//
//
//
AtlasMesh::CellsContainer::Pointer
MultiResolutionAtlasMesher
::GetCells( const AtlasMesh::PointsContainer* position ) const
{

#if 0
  // Write out a file that tetgen can read
  const std::string  tetgenInputFileName = "pointSet.node";
  std::ofstream  out( tetgenInputFileName.c_str() );
  if ( out.bad() )
  {
    itkExceptionMacro( << "Can't open " << tetgenInputFileName << " for writing." );
  }
  out << position->Size() << " 3  0  0" << std::endl;
  int counter = 0;
  for ( kvl::AtlasMesh::PointsContainer::ConstIterator  it = position->Begin();
        it != position->End(); ++it, ++counter )
  {
    out << counter << "  " << it.Value()[ 0 ] << "  " << it.Value()[ 1 ] << "  " << it.Value()[ 2 ] << std::endl;
  }
  out.close();
  std::cout << "Wrote " << tetgenInputFileName  << std::endl;

#if 0
  // Run tetgen
  std::cout << "Now do the following: \n" << std::endl;
  std::cout << "         ./tetgen " << tetgenInputFileName << std::endl;
  std::cout << "\nand then feed in any character to continue" << std::endl;
  char  dummy;
  std::cin >> dummy;
#else
  if ( execl( "tetgen", "tetgen", tetgenInputFileName.c_str(), 0 ) )
  {
    itkExceptionMacro( << "Couldn't start tetgen" );
  }
#endif


  // Read the output from tetgen and make a mesh from it
  const std::string  tetgenOutputFileName = "pointSet.1.ele";
  std::ifstream  in( tetgenOutputFileName.c_str() );
  if ( in.bad() )
  {
    itkExceptionMacro( << "Couldn't read from file " << tetgenOutputFileName );
  }

  typedef itk::AutomaticTopologyMeshSource< kvl::AtlasMesh >  MeshSourceType;
  MeshSourceType::Pointer  meshSource = MeshSourceType::New();
  const int size = 255;
  char buffer[ size ];
  in.getline( buffer, size ); // Skip first line
  while ( in.getline( buffer, size ) )
  {
    // Parse the line
    const std::string  line = buffer;
    std::istringstream  lineStream( line );
    unsigned long  tetrahedronId;
    unsigned long  point0Id;
    unsigned long  point1Id;
    unsigned long  point2Id;
    unsigned long  point3Id;
    lineStream >> tetrahedronId >> point0Id >> point1Id >> point2Id >> point3Id;

    meshSource->AddTetrahedron( point0Id,
                                point1Id,
                                point2Id,
                                point3Id );
  }
  in.close();


  return meshSource->GetOutput()->GetCells();

#else

#if 0
  // Write out a file that tetgen can read
  const std::string  tetgenInputFileName = "pointSet.node";
  std::ofstream  out( tetgenInputFileName.c_str() );
  if ( out.bad() )
  {
    itkExceptionMacro( << "Can't open " << tetgenInputFileName << " for writing." );
  }
  out << position->Size() << " 3  0  0" << std::endl;
  int counter = 0;
  for ( kvl::AtlasMesh::PointsContainer::ConstIterator  it = position->Begin();
        it != position->End(); ++it, ++counter )
  {
    out << counter << "  " << it.Value()[ 0 ] << "  " << it.Value()[ 1 ] << "  " << it.Value()[ 2 ] << std::endl;
  }
  out.close();
  std::cout << "Wrote " << tetgenInputFileName  << std::endl;

  // Run tetgen
  tetgenio  tetgenInput;
  if ( !tetgenInput.load_node( "pointSet" ) )
  {
    itkExceptionMacro( << "Can't open " << tetgenInputFileName << " for reading." );
  }
  //tetgenInput.save_nodes( "debugWithFile" );

#else

  // Initialize input structure for TetGen
  tetgenio  tetgenInput;
  tetgenInput.numberofpoints = position->Size();
  tetgenInput.mesh_dim = 3;
  tetgenInput.numberofpointattributes = 0;

  tetgenInput.pointlist = new REAL[ tetgenInput.numberofpoints * 3];
  if ( tetgenInput.pointlist == (REAL *) NULL )
  {
    itkExceptionMacro( << "Error:  Out of memory.\n" );
  }

  tetgenInput.firstnumber = 0;

  int counter = 0;
  for ( kvl::AtlasMesh::PointsContainer::ConstIterator  it = position->Begin();
        it != position->End(); ++it )
  {
    tetgenInput.pointlist[ counter++ ] = it.Value()[ 0 ];
    tetgenInput.pointlist[ counter++ ] = it.Value()[ 1 ];
    tetgenInput.pointlist[ counter++ ] = it.Value()[ 2 ];
  }

  //tetgenInput.save_nodes( "debugWithoutFile" );

#endif


  std::cout << "!!!!!!!!!!!!! Starting mesh generation " << std::endl;
  tetgenio  tetgenOutput;
  tetrahedralize( (char*)"", &tetgenInput, &tetgenOutput );
  std::cout << "!!!!!!!!!!!!! Finished mesh generation " << std::endl;


#if 0
  // Write tetgen output
  tetgenOutput.save_elements( "pointSet.1" );


  // Read the output from tetgen and make a mesh from it
  const std::string  tetgenOutputFileName = "pointSet.1.ele";
  std::ifstream  in( tetgenOutputFileName.c_str() );
  if ( in.bad() )
  {
    itkExceptionMacro( << "Couldn't read from file " << tetgenOutputFileName );
  }

  typedef itk::AutomaticTopologyMeshSource< kvl::AtlasMesh >  MeshSourceType;
  MeshSourceType::Pointer  meshSource = MeshSourceType::New();
  const int size = 255;
  char buffer[ size ];
  in.getline( buffer, size ); // Skip first line
  while ( in.getline( buffer, size ) )
  {
    // Parse the line
    const std::string  line = buffer;
    std::istringstream  lineStream( line );
    unsigned long  tetrahedronId;
    unsigned long  point0Id;
    unsigned long  point1Id;
    unsigned long  point2Id;
    unsigned long  point3Id;
    lineStream >> tetrahedronId >> point0Id >> point1Id >> point2Id >> point3Id;

    meshSource->AddTetrahedron( point0Id,
                                point1Id,
                                point2Id,
                                point3Id );
  }
  in.close();


  return meshSource->GetOutput()->GetCells();
#else

  MeshSourceType::Pointer  meshSource = MeshSourceType::New();
  for ( int i = 0; i < tetgenOutput.numberoftetrahedra; i++ )
  {
    unsigned long  point0Id = tetgenOutput.tetrahedronlist[ i * 4 + 0 ];
    unsigned long  point1Id = tetgenOutput.tetrahedronlist[ i * 4 + 1 ];
    unsigned long  point2Id = tetgenOutput.tetrahedronlist[ i * 4 + 2 ];
    unsigned long  point3Id = tetgenOutput.tetrahedronlist[ i * 4 + 3 ];

    meshSource->AddTetrahedron( point0Id,
                                point1Id,
                                point2Id,
                                point3Id );

#if 1
    {
      // Double-check that our tets are not negative volume.
      // Do this by calculating the volume of the tetrahedron; this should be positive.
      // In what follows, the matrix Lambda is the Jacobian of the transform from a standarized tetrahedron
      // ( ( 0 0 0 )^T, ( 1 0 0 )^T, ( 0 1 0 )^T, ( 0 0 1 )^T ), which has volume 1/6, to the actual tetrahedron
      const float x0 = position->ElementAt( point0Id )[ 0 ];
      const float y0 = position->ElementAt( point0Id )[ 1 ];
      const float z0 = position->ElementAt( point0Id )[ 2 ];

      const float x1 = position->ElementAt( point1Id )[ 0 ];
      const float y1 = position->ElementAt( point1Id )[ 1 ];
      const float z1 = position->ElementAt( point1Id )[ 2 ];

      const float x2 = position->ElementAt( point2Id )[ 0 ];
      const float y2 = position->ElementAt( point2Id )[ 1 ];
      const float z2 = position->ElementAt( point2Id )[ 2 ];

      const float x3 = position->ElementAt( point3Id )[ 0 ];
      const float y3 = position->ElementAt( point3Id )[ 1 ];
      const float z3 = position->ElementAt( point3Id )[ 2 ];

      const float  lambda11 = -x0 + x1;
      const float  lambda21 = -y0 + y1;
      const float  lambda31 = -z0 + z1;
      const float  lambda12 = -x0 + x2;
      const float  lambda22 = -y0 + y2;
      const float  lambda32 = -z0 + z2;
      const float  lambda13 = -x0 + x3;
      const float  lambda23 = -y0 + y3;
      const float  lambda33 = -z0 + z3;
      const float  volume = ( lambda11 * ( lambda22*lambda33 - lambda32*lambda23 )
                              - lambda12 * ( lambda21*lambda33 - lambda31*lambda23 )
                              + lambda13 * ( lambda21*lambda32 - lambda31*lambda22 ) ) / 6;
      if ( volume <= 0  )
      {
        std::cout << "****************************************" << std::endl;
        std::cout << "****************************************" << std::endl;
        std::cout << "Ouch: TetGen has generated a tetrahedron with negative volume! " << std::endl;
        std::cout << "      tetrahedronId: " << i << std::endl;
        std::cout << "                 p0: " << position->ElementAt( point0Id ) << std::endl;
        std::cout << "                 p1: " << position->ElementAt( point1Id ) << std::endl;
        std::cout << "                 p2: " << position->ElementAt( point2Id ) << std::endl;
        std::cout << "                 p3: " << position->ElementAt( point3Id ) << std::endl;
        std::cout << "****************************************" << std::endl;
        std::cout << "****************************************" << std::endl;

        exit( -1 );
      }

    }
#endif


  }

  return meshSource->GetOutput()->GetCells();

#endif


#endif

}



//
//
//
AtlasMeshCollection::Pointer
MultiResolutionAtlasMesher
::GetMeshCollection( AtlasMesh::PointsContainer* referencePosition,
                     std::vector< AtlasMesh::PointsContainer::Pointer >& positions,
                     float stiffness ) const
{
  // Construct the cells by running TetGen on the referencePosition point set
  AtlasMesh::CellsContainer::Pointer  cells = this->GetCells( referencePosition );


  // Also get point parameters. Assign flat alphas as a starting point. Vertices lying on the border
  // can not move freely and belong to
  // first class
  kvl::AtlasAlphasType   flatAlphasEntry( m_NumberOfClasses );
  flatAlphasEntry.Fill( 1.0f / static_cast< float >( m_NumberOfClasses ) );

  kvl::AtlasAlphasType   borderAlphasEntry( m_NumberOfClasses );
  borderAlphasEntry.Fill( 0.0f );
  borderAlphasEntry[ 0 ] = 1.0f;

  AtlasMesh::PointDataContainer::Pointer  pointData = AtlasMesh::PointDataContainer::New();
  for ( kvl::AtlasMesh::PointsContainer::ConstIterator  pointIt = referencePosition->Begin();
        pointIt != referencePosition->End();
        ++pointIt )
  {
    kvl::AtlasMesh::PixelType  pointParameters;

    pointParameters.m_Alphas = flatAlphasEntry;
    pointParameters.m_CanChangeAlphas = true;

    if ( ( fabs( pointIt.Value()[ 0 ] ) < 1e-3 ) || ( fabs( pointIt.Value()[ 0 ] - ( m_DomainSize[ 0 ] - 1 ) ) < 1e-3 ) )
    {
      pointParameters.m_CanMoveX = false;
      pointParameters.m_Alphas = borderAlphasEntry;
      pointParameters.m_CanChangeAlphas = false;
    }
    else
    {
      pointParameters.m_CanMoveX = true;
    }

    if ( ( fabs( pointIt.Value()[ 1 ] ) < 1e-3 ) || ( fabs( pointIt.Value()[ 1 ] - ( m_DomainSize[ 1 ] - 1 ) ) < 1e-3 ) )
    {
      pointParameters.m_CanMoveY = false;
      pointParameters.m_Alphas = borderAlphasEntry;
      pointParameters.m_CanChangeAlphas = false;
    }
    else
    {
      pointParameters.m_CanMoveY = true;
    }

    if ( ( fabs( pointIt.Value()[ 2 ] ) < 1e-3 ) || ( fabs( pointIt.Value()[ 2 ] - ( m_DomainSize[ 2 ] - 1 ) ) < 1e-3 ) )
    {
      pointParameters.m_CanMoveZ = false;
      pointParameters.m_Alphas = borderAlphasEntry;
      pointParameters.m_CanChangeAlphas = false;
    }
    else
    {
      pointParameters.m_CanMoveZ = true;
    }

#if 0
    {
      if ( pointIt.Index() == 166 )
      {
        std::cout << "@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@" << std::endl;
        std::cout << "Having point with index: " << pointIt.Index() << std::endl;
        std::cout <<  "           position is: " << pointIt.Value() << std::endl;
        std::cout <<  "           test for component [2]:" << ( fabs( pointIt.Value()[ 2 ] - ( m_DomainSize[ 2 ] - 1 ) ) < 1e-5 ) << std::endl;
        std::cout <<  "           fabs( pointIt.Value()[ 2 ] - ( m_DomainSize[ 2 ] - 1 ) ): "
                  << fabs( pointIt.Value()[ 2 ] - ( m_DomainSize[ 2 ] - 1 ) ) << std::endl;
        std::cout << "        pointParameters.m_CanMoveZ: " << pointParameters.m_CanMoveZ << std::endl;
        std::cout << "        pointIt.Value()[ 2 ]: " <<  pointIt.Value()[ 2 ] << std::endl;
        std::cout << "        ( m_DomainSize[ 2 ] - 1 ): " << ( m_DomainSize[ 2 ] - 1 ) << std::endl;
        std::cout << "@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@" << std::endl;
      }
    }
#endif

    pointData->InsertElement( pointIt.Index(), pointParameters );
  }


  // Now combine everything into a mesh collection
  AtlasMeshCollection::Pointer  meshCollection = AtlasMeshCollection::New();
  meshCollection->SetCells( cells );
  meshCollection->SetReferencePosition( referencePosition );
  meshCollection->SetPositions( positions );
  meshCollection->SetPointParameters( pointData );
  meshCollection->SetK( stiffness );


#if 1
  {
    //meshCollection->Write( "debug.txt" );
    //std::cout << "Enter any character to continue" << std::endl;
    //char dummy;
    //std::cin >> dummy;


    // Double-check that our tets are not negative volume.
    // Do this by calculating the volume of the tetrahedron; this should be positive.
    for ( int meshNumber = this->GetNumberOfMeshes(); meshNumber >= 0; meshNumber-- )
    {

      // Retrieve the corresponding original point set
      AtlasMesh::PointsContainer::Pointer  thisPosition;
      if ( meshNumber < static_cast< int >( this->GetNumberOfMeshes() ) )
      {
        thisPosition = meshCollection->GetPositions()[ meshNumber ];
      }
      else
      {
        thisPosition = meshCollection->GetReferencePosition();
      }


      // Loop over all tetrahedra
      for ( AtlasMesh::CellsContainer::ConstIterator  cellIt = meshCollection->GetCells()->Begin();
            cellIt != meshCollection->GetCells()->End(); ++cellIt )
      {
        const AtlasMesh::CellType*  cell = cellIt.Value();

        if ( cell->GetType() != AtlasMesh::CellType::TETRAHEDRON_CELL )
        {
          continue;
        }

        AtlasMesh::CellType::PointIdConstIterator  pit = cell->PointIdsBegin();
        AtlasMesh::CellIdentifier  point0Id = *pit;
        ++pit;
        AtlasMesh::CellIdentifier  point1Id = *pit;
        ++pit;
        AtlasMesh::CellIdentifier  point2Id = *pit;
        ++pit;
        AtlasMesh::CellIdentifier  point3Id = *pit;


        // In what follows, the matrix Lambda is the Jacobian of the transform from a standarized tetrahedron
        // ( ( 0 0 0 )^T, ( 1 0 0 )^T, ( 0 1 0 )^T, ( 0 0 1 )^T ), which has volume 1/6, to the actual tetrahedron
        const float x0 = thisPosition->ElementAt( point0Id )[ 0 ];
        const float y0 = thisPosition->ElementAt( point0Id )[ 1 ];
        const float z0 = thisPosition->ElementAt( point0Id )[ 2 ];

        const float x1 = thisPosition->ElementAt( point1Id )[ 0 ];
        const float y1 = thisPosition->ElementAt( point1Id )[ 1 ];
        const float z1 = thisPosition->ElementAt( point1Id )[ 2 ];

        const float x2 = thisPosition->ElementAt( point2Id )[ 0 ];
        const float y2 = thisPosition->ElementAt( point2Id )[ 1 ];
        const float z2 = thisPosition->ElementAt( point2Id )[ 2 ];

        const float x3 = thisPosition->ElementAt( point3Id )[ 0 ];
        const float y3 = thisPosition->ElementAt( point3Id )[ 1 ];
        const float z3 = thisPosition->ElementAt( point3Id )[ 2 ];

        const float  lambda11 = -x0 + x1;
        const float  lambda21 = -y0 + y1;
        const float  lambda31 = -z0 + z1;
        const float  lambda12 = -x0 + x2;
        const float  lambda22 = -y0 + y2;
        const float  lambda32 = -z0 + z2;
        const float  lambda13 = -x0 + x3;
        const float  lambda23 = -y0 + y3;
        const float  lambda33 = -z0 + z3;
        const float  volume = ( lambda11 * ( lambda22*lambda33 - lambda32*lambda23 )
                                - lambda12 * ( lambda21*lambda33 - lambda31*lambda23 )
                                + lambda13 * ( lambda21*lambda32 - lambda31*lambda22 ) ) / 6;
        if ( volume <= 0 )
        {
          std::cout << "****************************************" << std::endl;
          std::cout << "****************************************" << std::endl;
          std::cout << "Ouch: TetGen has generated a tetrahedron with negative volume! " << std::endl;
          if ( meshNumber < static_cast< int >( this->GetNumberOfMeshes() ) )
          {
            std::cout << "         meshNumber: " << meshNumber << std::endl;
          }
          else
          {
            std::cout << "      referenceMesh " << std::endl;
          }
          std::cout << "      tetrahedronId: " << cellIt.Index() << std::endl;
          std::cout << "               p0Id: " << point0Id << std::endl;
          std::cout << "               p1Id: " << point1Id << std::endl;
          std::cout << "               p2Id: " << point2Id << std::endl;
          std::cout << "               p3Id: " << point3Id << std::endl;
          std::cout << "                 p0: " << thisPosition->ElementAt( point0Id ) << std::endl;
          std::cout << "                 p1: " << thisPosition->ElementAt( point1Id ) << std::endl;
          std::cout << "                 p2: " << thisPosition->ElementAt( point2Id ) << std::endl;
          std::cout << "                 p3: " << thisPosition->ElementAt( point3Id ) << std::endl;
          std::cout << "             volume: " << volume << std::endl;
          std::cout << "****************************************" << std::endl;
          std::cout << "****************************************" << std::endl;

          exit( -1 );
        }

      } // End loop over all tetrahedra

    } // End loop over all meshes


  }
#endif




  return meshCollection;
}



//
//
//
void
MultiResolutionAtlasMesher
::Upsample()
{

  // Loop over all hexahedra
  std::vector< AtlasMesh::PointsContainer::Pointer >  upsampledPositions;
  for ( unsigned int i = 0; i < m_NumberOfMeshes; i++ )
  {
    upsampledPositions.push_back( AtlasMesh::PointsContainer::New() );
  }
  MeshSourceType::Pointer  upsampledReferencePositionSource = MeshSourceType::New();
  AtlasMesh::CellsContainer::Pointer  upsampledHexahedra = AtlasMesh::CellsContainer::New();

  float  precision[ 3 ];  // Used to convert the reference positions into exact integer coordinates
  for ( int i = 0; i < 3; i++ )
  {
    const int  factor = static_cast< int >( pow( 2, m_NumberOfUpsamplingSteps ) );
    const int  finalSize = factor * m_InitialSize[ i ] - ( factor - 1 );
    precision[ i ] = static_cast< float >( finalSize - 1 ) / static_cast< float >( m_DomainSize[ i ] - 1 );
  }

  for ( AtlasMesh::CellsContainer::ConstIterator  hexIt = m_Hexahedra->Begin();
        hexIt != m_Hexahedra->End(); ++hexIt )
  {
    // Retrieve the ids of the 8 corners
    const AtlasMesh::CellType*  cell = hexIt.Value();

    AtlasMesh::CellType::PointIdConstIterator  pit = cell->PointIdsBegin();
    const AtlasMesh::PointIdentifier  p0Id = *pit;
    ++pit;
    const AtlasMesh::PointIdentifier  p1Id = *pit;
    ++pit;
    const AtlasMesh::PointIdentifier  p2Id = *pit;
    ++pit;
    const AtlasMesh::PointIdentifier  p3Id = *pit;
    ++pit;
    const AtlasMesh::PointIdentifier  p4Id = *pit;
    ++pit;
    const AtlasMesh::PointIdentifier  p5Id = *pit;
    ++pit;
    const AtlasMesh::PointIdentifier  p6Id = *pit;
    ++pit;
    const AtlasMesh::PointIdentifier  p7Id = *pit;
    ++pit;


    // Clean up memory of this hexahedral cell
    delete hexIt.Value();


    // Decide whether or not we're gonna subdivide this hexahedron
    bool  subdivideHexahedron = true;
    if ( m_TryToBeSparse )
    {
      // Look up the label with the highest alpha in the first corner point
      int  maximumAlphaLabelNumber = 0;
      float  maximumAlpha = itk::NumericTraits< float >::min();
      for ( unsigned int classNumber = 0; classNumber < m_NumberOfClasses; classNumber++ )
      {
        if ( m_Current->GetPointParameters()->ElementAt( p0Id ).m_Alphas[ classNumber ] > maximumAlpha )
        {
          maximumAlpha = m_Current->GetPointParameters()->ElementAt( p0Id ).m_Alphas[ classNumber ];
          maximumAlphaLabelNumber = classNumber;
        }
      }


      // Look at the alphas in each of the corners points
      const float threshold = 0.90f;
      if ( ( m_Current->GetPointParameters()->ElementAt( p0Id ).m_Alphas[ maximumAlphaLabelNumber ] >= threshold ) &&
           ( m_Current->GetPointParameters()->ElementAt( p1Id ).m_Alphas[ maximumAlphaLabelNumber ] >= threshold ) &&
           ( m_Current->GetPointParameters()->ElementAt( p2Id ).m_Alphas[ maximumAlphaLabelNumber ] >= threshold ) &&
           ( m_Current->GetPointParameters()->ElementAt( p3Id ).m_Alphas[ maximumAlphaLabelNumber ] >= threshold ) &&
           ( m_Current->GetPointParameters()->ElementAt( p4Id ).m_Alphas[ maximumAlphaLabelNumber ] >= threshold ) &&
           ( m_Current->GetPointParameters()->ElementAt( p5Id ).m_Alphas[ maximumAlphaLabelNumber ] >= threshold ) &&
           ( m_Current->GetPointParameters()->ElementAt( p6Id ).m_Alphas[ maximumAlphaLabelNumber ] >= threshold ) &&
           ( m_Current->GetPointParameters()->ElementAt( p7Id ).m_Alphas[ maximumAlphaLabelNumber ] >= threshold ) )
      {
        subdivideHexahedron = false;
      }

    }



    if ( subdivideHexahedron )
    {

      // The values of these ids will be filled on from the reference mesh, and used by all other meshes
      AtlasMesh::PointIdentifier  p0UpsampledId;
      AtlasMesh::PointIdentifier  p1UpsampledId;
      AtlasMesh::PointIdentifier  p2UpsampledId;
      AtlasMesh::PointIdentifier  p3UpsampledId;
      AtlasMesh::PointIdentifier  p4UpsampledId;
      AtlasMesh::PointIdentifier  p5UpsampledId;
      AtlasMesh::PointIdentifier  p6UpsampledId;
      AtlasMesh::PointIdentifier  p7UpsampledId;
      AtlasMesh::PointIdentifier  p01UpsampledId;
      AtlasMesh::PointIdentifier  p02UpsampledId;
      AtlasMesh::PointIdentifier  p13UpsampledId;
      AtlasMesh::PointIdentifier  p23UpsampledId;
      AtlasMesh::PointIdentifier  p15UpsampledId;
      AtlasMesh::PointIdentifier  p37UpsampledId;
      AtlasMesh::PointIdentifier  p26UpsampledId;
      AtlasMesh::PointIdentifier  p04UpsampledId;
      AtlasMesh::PointIdentifier  p57UpsampledId;
      AtlasMesh::PointIdentifier  p67UpsampledId;
      AtlasMesh::PointIdentifier  p46UpsampledId;
      AtlasMesh::PointIdentifier  p45UpsampledId;
      AtlasMesh::PointIdentifier  p0123UpsampledId;
      AtlasMesh::PointIdentifier  p1357UpsampledId;
      AtlasMesh::PointIdentifier  p2367UpsampledId;
      AtlasMesh::PointIdentifier  p0246UpsampledId;
      AtlasMesh::PointIdentifier  p0145UpsampledId;
      AtlasMesh::PointIdentifier  p4567UpsampledId;
      AtlasMesh::PointIdentifier  pMiddleUpsampledId;

      // Loop over all meshes, starting with the reference mesh to obtain the correct ids
      for ( int meshNumber=this->GetNumberOfMeshes(); meshNumber >= 0; meshNumber-- )
      {

        // Retrieve the corresponding original point set
        AtlasMesh::PointsContainer::Pointer  thisPosition;
        if ( meshNumber < static_cast< int >( this->GetNumberOfMeshes() ) )
        {
          thisPosition = m_Current->GetPositions()[ meshNumber ];
        }
        else
        {
          thisPosition = m_Current->GetReferencePosition();
        }


        // Retrieve the position of the 8 corners
        AtlasMesh::PointType  p0 = thisPosition->ElementAt( p0Id );
        AtlasMesh::PointType  p1 = thisPosition->ElementAt( p1Id );
        AtlasMesh::PointType  p2 = thisPosition->ElementAt( p2Id );
        AtlasMesh::PointType  p3 = thisPosition->ElementAt( p3Id );
        AtlasMesh::PointType  p4 = thisPosition->ElementAt( p4Id );
        AtlasMesh::PointType  p5 = thisPosition->ElementAt( p5Id );
        AtlasMesh::PointType  p6 = thisPosition->ElementAt( p6Id );
        AtlasMesh::PointType  p7 = thisPosition->ElementAt( p7Id );


        // Construct points on the middle of each cube edge
        AtlasMesh::PointType  p01;
        AtlasMesh::PointType  p02;
        AtlasMesh::PointType  p13;
        AtlasMesh::PointType  p23;
        AtlasMesh::PointType  p15;
        AtlasMesh::PointType  p37;
        AtlasMesh::PointType  p26;
        AtlasMesh::PointType  p04;
        AtlasMesh::PointType  p57;
        AtlasMesh::PointType  p67;
        AtlasMesh::PointType  p46;
        AtlasMesh::PointType  p45;
        AtlasMesh::PointType  p0123;
        AtlasMesh::PointType  p1357;
        AtlasMesh::PointType  p2367;
        AtlasMesh::PointType  p0246;
        AtlasMesh::PointType  p0145;
        AtlasMesh::PointType  p4567;
        AtlasMesh::PointType  pMiddle;
        if ( meshNumber == static_cast< int >( this->GetNumberOfMeshes() ) )
        {
          // Convert the input into integers
          for ( int i = 0; i < 3; i++ )
          {
            p0[ i ] = static_cast< int >( p0[ i ] * precision[ i ] + 0.5 );
            p1[ i ] = static_cast< int >( p1[ i ] * precision[ i ] + 0.5 );
            p2[ i ] = static_cast< int >( p2[ i ] * precision[ i ] + 0.5 );
            p3[ i ] = static_cast< int >( p3[ i ] * precision[ i ] + 0.5 );
            p4[ i ] = static_cast< int >( p4[ i ] * precision[ i ] + 0.5 );
            p5[ i ] = static_cast< int >( p5[ i ] * precision[ i ] + 0.5 );
            p6[ i ] = static_cast< int >( p6[ i ] * precision[ i ] + 0.5 );
            p7[ i ] = static_cast< int >( p7[ i ] * precision[ i ] + 0.5 );
          }

          this->GetUpsampledHexahedronPoints< int >( p0, p1, p2, p3, p4, p5, p6, p7,
              p01, p02, p13, p23, p15, p37, p26, p04, p57, p67, p46, p45,
              p0123, p1357, p2367, p0246, p0145, p4567,
              pMiddle );
        }
        else
        {
          this->GetUpsampledHexahedronPoints< float >( p0, p1, p2, p3, p4, p5, p6, p7,
              p01, p02, p13, p23, p15, p37, p26, p04, p57, p67, p46, p45,
              p0123, p1357, p2367, p0246, p0145, p4567,
              pMiddle );
        }


        // If this is the reference mesh, add the points while simulatenously creating the correct hexahedra
        // and looking up the ids of the points to be used for the other meshes. Otherwise, just add the points.
        if ( meshNumber == static_cast< int >( this->GetNumberOfMeshes() ) )
        {
          // Create 8 sub-hexahedra
          this->AddHexahedron( upsampledReferencePositionSource,
                               upsampledHexahedra,
                               p0, p01, p02, p0123,
                               p04, p0145, p0246, pMiddle,
                               p0UpsampledId, p01UpsampledId, p02UpsampledId, p0123UpsampledId,
                               p04UpsampledId, p0145UpsampledId, p0246UpsampledId, pMiddleUpsampledId );
          this->AddHexahedron( upsampledReferencePositionSource,
                               upsampledHexahedra,
                               p01, p1, p0123, p13,
                               p0145, p15, pMiddle, p1357,
                               p01UpsampledId, p1UpsampledId, p0123UpsampledId, p13UpsampledId,
                               p0145UpsampledId, p15UpsampledId, pMiddleUpsampledId, p1357UpsampledId );
          this->AddHexahedron( upsampledReferencePositionSource,
                               upsampledHexahedra,
                               p02, p0123, p2, p23,
                               p0246, pMiddle, p26, p2367,
                               p02UpsampledId, p0123UpsampledId, p2UpsampledId, p23UpsampledId,
                               p0246UpsampledId, pMiddleUpsampledId, p26UpsampledId, p2367UpsampledId );
          this->AddHexahedron( upsampledReferencePositionSource,
                               upsampledHexahedra,
                               p0123, p13, p23, p3,
                               pMiddle, p1357, p2367, p37,
                               p0123UpsampledId, p13UpsampledId, p23UpsampledId, p3UpsampledId,
                               pMiddleUpsampledId, p1357UpsampledId, p2367UpsampledId, p37UpsampledId );
          this->AddHexahedron( upsampledReferencePositionSource,
                               upsampledHexahedra,
                               p04, p0145, p0246, pMiddle,
                               p4, p45, p46, p4567,
                               p04UpsampledId, p0145UpsampledId, p0246UpsampledId, pMiddleUpsampledId,
                               p4UpsampledId, p45UpsampledId, p46UpsampledId, p4567UpsampledId );
          this->AddHexahedron( upsampledReferencePositionSource,
                               upsampledHexahedra,
                               p0145, p15, pMiddle, p1357,
                               p45, p5, p4567, p57,
                               p0145UpsampledId, p15UpsampledId, pMiddleUpsampledId, p1357UpsampledId,
                               p45UpsampledId, p5UpsampledId, p4567UpsampledId, p57UpsampledId );
          this->AddHexahedron( upsampledReferencePositionSource,
                               upsampledHexahedra,
                               p0246, pMiddle, p26, p2367,
                               p46, p4567, p6, p67,
                               p0246UpsampledId, pMiddleUpsampledId, p26UpsampledId, p2367UpsampledId,
                               p46UpsampledId, p4567UpsampledId, p6UpsampledId, p67UpsampledId );
          this->AddHexahedron( upsampledReferencePositionSource,
                               upsampledHexahedra,
                               pMiddle, p1357, p2367, p37,
                               p4567, p57, p67, p7,
                               pMiddleUpsampledId, p1357UpsampledId, p2367UpsampledId, p37UpsampledId,
                               p4567UpsampledId, p57UpsampledId, p67UpsampledId, p7UpsampledId );
        }
        else
        {
          upsampledPositions[ meshNumber ]->InsertElement( p0UpsampledId, p0 );
          upsampledPositions[ meshNumber ]->InsertElement( p1UpsampledId, p1 );
          upsampledPositions[ meshNumber ]->InsertElement( p2UpsampledId, p2 );
          upsampledPositions[ meshNumber ]->InsertElement( p3UpsampledId, p3 );
          upsampledPositions[ meshNumber ]->InsertElement( p4UpsampledId, p4 );
          upsampledPositions[ meshNumber ]->InsertElement( p5UpsampledId, p5 );
          upsampledPositions[ meshNumber ]->InsertElement( p6UpsampledId, p6 );
          upsampledPositions[ meshNumber ]->InsertElement( p7UpsampledId, p7 );

          upsampledPositions[ meshNumber ]->InsertElement( p01UpsampledId, p01 );
          upsampledPositions[ meshNumber ]->InsertElement( p02UpsampledId, p02 );
          upsampledPositions[ meshNumber ]->InsertElement( p13UpsampledId, p13 );
          upsampledPositions[ meshNumber ]->InsertElement( p23UpsampledId, p23 );
          upsampledPositions[ meshNumber ]->InsertElement( p15UpsampledId, p15 );
          upsampledPositions[ meshNumber ]->InsertElement( p37UpsampledId, p37 );
          upsampledPositions[ meshNumber ]->InsertElement( p26UpsampledId, p26 );
          upsampledPositions[ meshNumber ]->InsertElement( p04UpsampledId, p04 );
          upsampledPositions[ meshNumber ]->InsertElement( p57UpsampledId, p57 );
          upsampledPositions[ meshNumber ]->InsertElement( p67UpsampledId, p67 );
          upsampledPositions[ meshNumber ]->InsertElement( p46UpsampledId, p46 );
          upsampledPositions[ meshNumber ]->InsertElement( p45UpsampledId, p45 );

          upsampledPositions[ meshNumber ]->InsertElement( p0123UpsampledId, p0123 );
          upsampledPositions[ meshNumber ]->InsertElement( p1357UpsampledId, p1357 );
          upsampledPositions[ meshNumber ]->InsertElement( p2367UpsampledId, p2367 );
          upsampledPositions[ meshNumber ]->InsertElement( p0246UpsampledId, p0246 );
          upsampledPositions[ meshNumber ]->InsertElement( p0145UpsampledId, p0145 );
          upsampledPositions[ meshNumber ]->InsertElement( p4567UpsampledId, p4567 );

          upsampledPositions[ meshNumber ]->InsertElement( pMiddleUpsampledId, pMiddle );
        }

      } // End loop over all meshes


    }
    else
    {
      // Don't subdivide the hexahedron; just make it anew

      // The values of these ids will be filled on from the reference mesh, and used by all other meshes
      AtlasMesh::PointIdentifier  p0UpsampledId;
      AtlasMesh::PointIdentifier  p1UpsampledId;
      AtlasMesh::PointIdentifier  p2UpsampledId;
      AtlasMesh::PointIdentifier  p3UpsampledId;
      AtlasMesh::PointIdentifier  p4UpsampledId;
      AtlasMesh::PointIdentifier  p5UpsampledId;
      AtlasMesh::PointIdentifier  p6UpsampledId;
      AtlasMesh::PointIdentifier  p7UpsampledId;

      // Loop over all meshes
      for ( int meshNumber=this->GetNumberOfMeshes(); meshNumber >= 0; meshNumber-- )
      {

        // Retrieve the corresponding original point set
        AtlasMesh::PointsContainer::Pointer  thisPosition;
        if ( meshNumber < static_cast< int >( this->GetNumberOfMeshes() ) )
        {
          thisPosition = m_Current->GetPositions()[ meshNumber ];
        }
        else
        {
          thisPosition = m_Current->GetReferencePosition();
        }


        // Retrieve the position of the 8 corners
        AtlasMesh::PointType  p0 = thisPosition->ElementAt( p0Id );
        AtlasMesh::PointType  p1 = thisPosition->ElementAt( p1Id );
        AtlasMesh::PointType  p2 = thisPosition->ElementAt( p2Id );
        AtlasMesh::PointType  p3 = thisPosition->ElementAt( p3Id );
        AtlasMesh::PointType  p4 = thisPosition->ElementAt( p4Id );
        AtlasMesh::PointType  p5 = thisPosition->ElementAt( p5Id );
        AtlasMesh::PointType  p6 = thisPosition->ElementAt( p6Id );
        AtlasMesh::PointType  p7 = thisPosition->ElementAt( p7Id );

        if ( meshNumber == static_cast< int >( this->GetNumberOfMeshes() ) )
        {
          // Convert the input into integers
          for ( int i = 0; i < 3; i++ )
          {
            p0[ i ] = static_cast< int >( p0[ i ] * precision[ i ] + 0.5 );
            p1[ i ] = static_cast< int >( p1[ i ] * precision[ i ] + 0.5 );
            p2[ i ] = static_cast< int >( p2[ i ] * precision[ i ] + 0.5 );
            p3[ i ] = static_cast< int >( p3[ i ] * precision[ i ] + 0.5 );
            p4[ i ] = static_cast< int >( p4[ i ] * precision[ i ] + 0.5 );
            p5[ i ] = static_cast< int >( p5[ i ] * precision[ i ] + 0.5 );
            p6[ i ] = static_cast< int >( p6[ i ] * precision[ i ] + 0.5 );
            p7[ i ] = static_cast< int >( p7[ i ] * precision[ i ] + 0.5 );
          }

          this->AddHexahedron( upsampledReferencePositionSource, upsampledHexahedra,
                               p0, p1, p2, p3, p4, p5, p6, p7,
                               p0UpsampledId, p1UpsampledId, p2UpsampledId, p3UpsampledId,
                               p4UpsampledId, p5UpsampledId, p6UpsampledId, p7UpsampledId );
        }
        else
        {
          upsampledPositions[ meshNumber ]->InsertElement( p0UpsampledId, p0 );
          upsampledPositions[ meshNumber ]->InsertElement( p1UpsampledId, p1 );
          upsampledPositions[ meshNumber ]->InsertElement( p2UpsampledId, p2 );
          upsampledPositions[ meshNumber ]->InsertElement( p3UpsampledId, p3 );
          upsampledPositions[ meshNumber ]->InsertElement( p4UpsampledId, p4 );
          upsampledPositions[ meshNumber ]->InsertElement( p5UpsampledId, p5 );
          upsampledPositions[ meshNumber ]->InsertElement( p6UpsampledId, p6 );
          upsampledPositions[ meshNumber ]->InsertElement( p7UpsampledId, p7 );

        }


      } // End loop over all meshes


    } // End test whether or not to split this tetrahedron

  } // End loop over all hexahedra


  // Retrieve the upsampled reference position
  AtlasMesh::PointsContainer::Pointer  upsampledReferencePosition
  = upsampledReferencePositionSource->GetOutput()->GetPoints();
  for ( AtlasMesh::PointsContainer::Iterator  it = upsampledReferencePosition->Begin();
        it != upsampledReferencePosition->End(); ++it )
  {
    for ( int i = 0; i < 3; i++ )
    {
      it.Value()[ i ] /= precision[ i ];
    }
  }


  // OK, so now we have a new reference position and position. Let's generate a mesh collection from that.
  m_Current = this->GetMeshCollection( upsampledReferencePosition, upsampledPositions, m_Current->GetK() );

  // Now initialize the estimator with the mesh
  m_Estimator->SetInitialMeshCollection( m_Current );

  // Also replace the old m_Hexahedra with the upsampeld one
  m_Hexahedra = upsampledHexahedra;

}




//
//
//
void
MultiResolutionAtlasMesher
::AddHexahedron( MeshSourceType* meshSource,
                 AtlasMesh::CellsContainer* hexahedra,
                 const AtlasMesh::PointType&  p0,
                 const AtlasMesh::PointType&  p1,
                 const AtlasMesh::PointType&  p2,
                 const AtlasMesh::PointType&  p3,
                 const AtlasMesh::PointType&  p4,
                 const AtlasMesh::PointType&  p5,
                 const AtlasMesh::PointType&  p6,
                 const AtlasMesh::PointType&  p7,
                 AtlasMesh::PointIdentifier&  p0Id,
                 AtlasMesh::PointIdentifier&  p1Id,
                 AtlasMesh::PointIdentifier&  p2Id,
                 AtlasMesh::PointIdentifier&  p3Id,
                 AtlasMesh::PointIdentifier&  p4Id,
                 AtlasMesh::PointIdentifier&  p5Id,
                 AtlasMesh::PointIdentifier&  p6Id,
                 AtlasMesh::PointIdentifier&  p7Id )
{

  // Retrieve the ids of these eight corner points. If they're new, a new id will automatically
  // be created
  p0Id = meshSource->AddPoint( p0 );
  p1Id = meshSource->AddPoint( p1 );
  p2Id = meshSource->AddPoint( p2 );
  p3Id = meshSource->AddPoint( p3 );
  p4Id = meshSource->AddPoint( p4 );
  p5Id = meshSource->AddPoint( p5 );
  p6Id = meshSource->AddPoint( p6 );
  p7Id = meshSource->AddPoint( p7 );


  if ( hexahedra )
  {
    // Create a hexahedral element for later usage
    typedef itk::HexahedronCell< AtlasMesh::CellType >  HexahedronCell;
    AtlasMesh::CellAutoPointer newCell;
    newCell.TakeOwnership( new HexahedronCell );
    newCell->SetPointId( 0, p0Id );
    newCell->SetPointId( 1, p1Id );
    newCell->SetPointId( 2, p2Id );
    newCell->SetPointId( 3, p3Id );
    newCell->SetPointId( 4, p4Id );
    newCell->SetPointId( 5, p5Id );
    newCell->SetPointId( 6, p6Id );
    newCell->SetPointId( 7, p7Id );

    // Add the cell
    hexahedra->InsertElement( hexahedra->Size(), newCell.ReleaseOwnership() );
  }


}



//
//
//
template<class TCoordRep>
void
MultiResolutionAtlasMesher
::GetUpsampledHexahedronPoints( const AtlasMesh::PointType&  p0,
                                const AtlasMesh::PointType&  p1,
                                const AtlasMesh::PointType&  p2,
                                const AtlasMesh::PointType&  p3,
                                const AtlasMesh::PointType&  p4,
                                const AtlasMesh::PointType&  p5,
                                const AtlasMesh::PointType&  p6,
                                const AtlasMesh::PointType&  p7,
                                AtlasMesh::PointType&  p01,
                                AtlasMesh::PointType&  p02,
                                AtlasMesh::PointType&  p13,
                                AtlasMesh::PointType&  p23,
                                AtlasMesh::PointType&  p15,
                                AtlasMesh::PointType&  p37,
                                AtlasMesh::PointType&  p26,
                                AtlasMesh::PointType&  p04,
                                AtlasMesh::PointType&  p57,
                                AtlasMesh::PointType&  p67,
                                AtlasMesh::PointType&  p46,
                                AtlasMesh::PointType&  p45,
                                AtlasMesh::PointType&  p0123,
                                AtlasMesh::PointType&  p1357,
                                AtlasMesh::PointType&  p2367,
                                AtlasMesh::PointType&  p0246,
                                AtlasMesh::PointType&  p0145,
                                AtlasMesh::PointType&  p4567,
                                AtlasMesh::PointType&  pMiddle )
{
  // Floating-point operations somehow don't seem to be entirely reproducible: we'd of course
  // like middle points calculated on adjacent hexahedra to be *exactly* the same, but using
  // normal floating-point ops, there sometimes appears to be a tiny difference in the calculated
  // location, resulting in the same point being split into two extremely close ones. To avoid this,
  // let's first convert to integer values, do the calculations on integers, and then convert back

  typedef itk::Point< TCoordRep, 3 >  InternalPointType;
  InternalPointType  internalP0;
  InternalPointType  internalP1;
  InternalPointType  internalP2;
  InternalPointType  internalP3;
  InternalPointType  internalP4;
  InternalPointType  internalP5;
  InternalPointType  internalP6;
  InternalPointType  internalP7;
  for ( int i = 0; i < 3; i++ )
  {
    internalP0[ i ] = static_cast< TCoordRep >( p0[ i ] );
    internalP1[ i ] = static_cast< TCoordRep >( p1[ i ] );
    internalP2[ i ] = static_cast< TCoordRep >( p2[ i ] );
    internalP3[ i ] = static_cast< TCoordRep >( p3[ i ] );
    internalP4[ i ] = static_cast< TCoordRep >( p4[ i ] );
    internalP5[ i ] = static_cast< TCoordRep >( p5[ i ] );
    internalP6[ i ] = static_cast< TCoordRep >( p6[ i ] );
    internalP7[ i ] = static_cast< TCoordRep >( p7[ i ] );
  }

  // std::cout << "internalP0: " << internalP0 << std::endl;
  // std::cout << "internalP1: " << internalP1 << std::endl;
  // std::cout << "internalP2: " << internalP2 << std::endl;
  // std::cout << "internalP3: " << internalP3 << std::endl;
  // std::cout << "internalP4: " << internalP4 << std::endl;
  // std::cout << "internalP5: " << internalP5 << std::endl;
  // std::cout << "internalP6: " << internalP6 << std::endl;
  // std::cout << "internalP7: " << internalP7 << std::endl;




  // Now calculate the position of the additional points
  for ( int i = 0; i < 3; i++ )
  {
#if 0
    // Points on the middle of each cube edge
    p01[ i ] = ( ( internalP0[ i ] + internalP1[ i ] ) / 2 ) / precision[ i ];
    p02[ i ] = ( ( internalP0[ i ] + internalP2[ i ] ) / 2 ) / precision[ i ];
    p13[ i ] = ( ( internalP1[ i ] + internalP3[ i ] ) / 2 ) / precision[ i ];
    p23[ i ] = ( ( internalP2[ i ] + internalP3[ i ] ) / 2 ) / precision[ i ];
    p15[ i ] = ( ( internalP1[ i ] + internalP5[ i ] ) / 2 ) / precision[ i ];
    p37[ i ] = ( ( internalP3[ i ] + internalP7[ i ] ) / 2 ) / precision[ i ];
    p26[ i ] = ( ( internalP2[ i ] + internalP6[ i ] ) / 2 ) / precision[ i ];
    p04[ i ] = ( ( internalP0[ i ] + internalP4[ i ] ) / 2 ) / precision[ i ];
    p57[ i ] = ( ( internalP5[ i ] + internalP7[ i ] ) / 2 ) / precision[ i ];
    p67[ i ] = ( ( internalP6[ i ] + internalP7[ i ] ) / 2 ) / precision[ i ];
    p46[ i ] = ( ( internalP4[ i ] + internalP6[ i ] ) / 2 ) / precision[ i ];
    p45[ i ] = ( ( internalP4[ i ] + internalP5[ i ] ) / 2 ) / precision[ i ];

    // Points in the middle of each cube face
    p0123[ i ] = ( ( internalP0[ i ] + internalP1[ i ] + internalP2[ i ] + internalP3[ i ] ) / 4 ) / precision[ i ];
    p1357[ i ] = ( ( internalP1[ i ] + internalP3[ i ] + internalP5[ i ] + internalP7[ i ] ) / 4 ) / precision[ i ];
    p2367[ i ] = ( ( internalP2[ i ] + internalP3[ i ] + internalP6[ i ] + internalP7[ i ] ) / 4 ) / precision[ i ];
    p0246[ i ] = ( ( internalP0[ i ] + internalP2[ i ] + internalP4[ i ] + internalP6[ i ] ) / 4 ) / precision[ i ];
    p0145[ i ] = ( ( internalP0[ i ] + internalP1[ i ] + internalP4[ i ] + internalP5[ i ] ) / 4 ) / precision[ i ];
    p4567[ i ] = ( ( internalP4[ i ] + internalP5[ i ] + internalP6[ i ] + internalP7[ i ] ) / 4 ) / precision[ i ];

    // Point in the middle of the cube
    pMiddle[ i ] = ( ( internalP0[ i ] + internalP1[ i ] + internalP2[ i ] + internalP3[ i ] + internalP4[ i ] + internalP5[ i ] + internalP6[ i ] + internalP7[ i ] ) / 8 ) / precision[ i ];


#else

    // Points on the middle of each cube edge
    p01[ i ] = ( internalP0[ i ] + internalP1[ i ] ) / 2;
    p02[ i ] = ( internalP0[ i ] + internalP2[ i ] ) / 2;
    p13[ i ] = ( internalP1[ i ] + internalP3[ i ] ) / 2;
    p23[ i ] = ( internalP2[ i ] + internalP3[ i ] ) / 2;
    p15[ i ] = ( internalP1[ i ] + internalP5[ i ] ) / 2;
    p37[ i ] = ( internalP3[ i ] + internalP7[ i ] ) / 2;
    p26[ i ] = ( internalP2[ i ] + internalP6[ i ] ) / 2;
    p04[ i ] = ( internalP0[ i ] + internalP4[ i ] ) / 2;
    p57[ i ] = ( internalP5[ i ] + internalP7[ i ] ) / 2;
    p67[ i ] = ( internalP6[ i ] + internalP7[ i ] ) / 2;
    p46[ i ] = ( internalP4[ i ] + internalP6[ i ] ) / 2;
    p45[ i ] = ( internalP4[ i ] + internalP5[ i ] ) / 2;

    // Points in the middle of each cube face
    p0123[ i ] = ( internalP0[ i ] + internalP1[ i ] + internalP2[ i ] + internalP3[ i ] ) / 4;
    p1357[ i ] = ( internalP1[ i ] + internalP3[ i ] + internalP5[ i ] + internalP7[ i ] ) / 4;
    p2367[ i ] = ( internalP2[ i ] + internalP3[ i ] + internalP6[ i ] + internalP7[ i ] ) / 4;
    p0246[ i ] = ( internalP0[ i ] + internalP2[ i ] + internalP4[ i ] + internalP6[ i ] ) / 4;
    p0145[ i ] = ( internalP0[ i ] + internalP1[ i ] + internalP4[ i ] + internalP5[ i ] ) / 4;
    p4567[ i ] = ( internalP4[ i ] + internalP5[ i ] + internalP6[ i ] + internalP7[ i ] ) / 4;

    // Point in the middle of the cube
    pMiddle[ i ] = ( internalP0[ i ] + internalP1[ i ] + internalP2[ i ] + internalP3[ i ] + internalP4[ i ] + internalP5[ i ] + internalP6[ i ] + internalP7[ i ] ) / 8;
#endif

  }


  // std::cout << "p0: " << p0 << std::endl;
  // std::cout << "p1: " << p1 << std::endl;
  // std::cout << "p2: " << p2 << std::endl;
  // std::cout << "p3: " << p3 << std::endl;
  // std::cout << "p4: " << p4 << std::endl;
  // std::cout << "p5: " << p5 << std::endl;
  // std::cout << "p6: " << p6 << std::endl;
  // std::cout << "p7: " << p7 << std::endl;
  // std::cout << "p01: " << p01 << std::endl;
  // std::cout << "   : " <<  p02  << std::endl;
  // std::cout << "   : " <<  p13  << std::endl;
  // std::cout << "   : " <<  p23  << std::endl;
  // std::cout << "   : " <<  p15  << std::endl;
  // std::cout << "   : " <<  p37  << std::endl;
  // std::cout << "   : " <<  p26  << std::endl;
  // std::cout << "   : " <<  p04  << std::endl;
  // std::cout << "   : " <<  p57  << std::endl;
  // std::cout << "   : " <<  p67  << std::endl;
  // std::cout << "   : " <<  p46  << std::endl;
  // std::cout << "   : " <<  p45  << std::endl;
  // std::cout << "   : " <<  p0123  << std::endl;
  // std::cout << "   : " <<  p1357  << std::endl;
  // std::cout << "   : " <<  p2367  << std::endl;
  // std::cout << "   : " <<  p0246  << std::endl;
  // std::cout << "   : " <<  p0145  << std::endl;
  // std::cout << "   : " <<  p4567  << std::endl;
  // std::cout << "   : " <<  pMiddle  << std::endl;



}




} // end namespace kvl
