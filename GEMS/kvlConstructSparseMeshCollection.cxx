/**
 * @file  kvlConstructSparseMeshCollection.cxx
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
#include "kvlAtlasParameterEstimator.h"
#include "kvlCompressionLookupTable.h"
#include "itkImageFileReader.h"
#include <fstream>


int main( int argc, char** argv )
{

  // Sanity check on input
  if ( argc < 6 )
  {
    std::cerr << "Usage: " << argv[ 0 ] << " meshSizeX meshSizeY meshSizeZ stiffness fileName1 [ fileName2 ... ]" << std::endl;

    return -1;
  }

  // Retrieve the input parameters
  std::ostringstream  inputParserStream;
  for ( int argumentNumber = 1; argumentNumber < 5; argumentNumber++ )
  {
    inputParserStream << argv[ argumentNumber ] << " ";
  }
  std::istringstream  inputStream( inputParserStream.str().c_str() );
  unsigned int  meshSizeX;
  unsigned int  meshSizeY;
  unsigned int  meshSizeZ;
  float  stiffness;
  inputStream >> meshSizeX >> meshSizeY >> meshSizeZ >> stiffness;

  try
  {
    // Read the images (in original ushort format)
    typedef kvl::CompressionLookupTable::ImageType  InputImageType;
    std::vector< InputImageType::ConstPointer >  originalImages;
    for ( int argumentNumber = 5; argumentNumber < argc; argumentNumber++ )
    {
      // Read the input image
      typedef itk::ImageFileReader< InputImageType >  ReaderType;
      ReaderType::Pointer  reader = ReaderType::New();
      reader->SetFileName( argv[ argumentNumber ] );
      reader->Update();
      InputImageType::ConstPointer  originalImage = reader->GetOutput();

      // Over-ride the spacing and origin since at this point we can't deal with that
      const double spacing[] = { 1, 1, 1 };
      const double origin[] = { 0, 0, 0 };
      const_cast< InputImageType* >( originalImage.GetPointer() )->SetSpacing( spacing );
      const_cast< InputImageType* >( originalImage.GetPointer() )->SetOrigin( origin );

      // Remember this image
      originalImages.push_back( originalImage );
    }

    // Build a lookup table that maps the original intensities onto uchar starting
    // at 0 and densely packed
    kvl::CompressionLookupTable::Pointer  compressor = kvl::CompressionLookupTable::New();
    compressor->Construct( originalImages );
    compressor->Write( "compressionLookupTable.txt" );

    // Collect the label images resulting from pushing the original images through the
    // lookup table
    typedef kvl::CompressionLookupTable::CompressedImageType  OutputImageType;
    std::vector< OutputImageType::ConstPointer >  labelImages;
    for ( std::vector< InputImageType::ConstPointer >::const_iterator it = originalImages.begin();
          it != originalImages.end(); ++it )
    {
      labelImages.push_back( ( compressor->CompressImage( ( *it ).GetPointer() ) ).GetPointer() );
    }


    // Set up an estimator and a mesh collection to go with it
    kvl::AtlasParameterEstimator::Pointer  estimator = kvl::AtlasParameterEstimator::New();
    estimator->SetLabelImages( labelImages );

    unsigned int  meshSize[ 3 ];
    meshSize[ 0 ] = static_cast< unsigned int >( meshSizeX );
    meshSize[ 1 ] = static_cast< unsigned int >( meshSizeY );
    meshSize[ 2 ] = static_cast< unsigned int >( meshSizeZ );

    kvl::AtlasParameterEstimator::LabelImageType::SizeType  labelImageSize =
      estimator->GetLabelImage( 0 )->GetLargestPossibleRegion().GetSize();
    unsigned int  domainSize[ 3 ];
    domainSize[ 0 ] = static_cast< unsigned int >( labelImageSize[ 0 ] );
    domainSize[ 1 ] = static_cast< unsigned int >( labelImageSize[ 1 ] );
    domainSize[ 2 ] = static_cast< unsigned int >( labelImageSize[ 2 ] );

    const unsigned int  numberOfClasses = estimator->GetNumberOfClasses();

    const unsigned int  numberOfMeshes = estimator->GetNumberOfLabelImages();

    kvl::AtlasMeshCollection::Pointer  meshCollection = kvl::AtlasMeshCollection::New();
    meshCollection->Construct( meshSize, domainSize, stiffness,
                               numberOfClasses, numberOfMeshes );

    // Estimate the mesh with infinite stiffness
    meshCollection->SetK( 1000 );
    estimator->SetInitialMeshCollection( meshCollection );
    estimator->Estimate();
    meshCollection->Write( "original.txt" );


    // Now collect a set of points that do NOT lie on an edge with different alphas in its edge points
    kvl::AtlasMesh::PointsContainer::Pointer  sparsePoints = kvl::AtlasMesh::PointsContainer::New();
    for ( kvl::AtlasMesh::PointsContainer::ConstIterator  it = meshCollection->GetReferencePosition()->Begin();
          it != meshCollection->GetReferencePosition()->End(); ++it )
    {
      sparsePoints->InsertElement( it.Index(), it.Value() );
    }

    for ( kvl::AtlasMesh::CellsContainer::ConstIterator  cellIt = meshCollection->GetCells()->Begin();
          cellIt != meshCollection->GetCells()->End(); ++cellIt )
    {
      const kvl::AtlasMesh::CellType*  cell = cellIt.Value();

      if ( cell->GetType() != kvl::AtlasMesh::CellType::LINE_CELL )
      {
        continue;
      }

      // Get the point ids of the two edge points
      kvl::AtlasMesh::CellType::PointIdConstIterator  pointIt = cell->PointIdsBegin();
      const kvl::AtlasMesh::PointIdentifier  point0Id = *pointIt;
      ++pointIt;
      const kvl::AtlasMesh::PointIdentifier  point1Id = *pointIt;

      // Check the alphas in the two points
      const float  threshold = 0.9999f;
      if ( ( meshCollection->GetPointParameters()->ElementAt( point0Id ).m_Alphas[ 0 ] >= threshold ) &&
           ( meshCollection->GetPointParameters()->ElementAt( point1Id ).m_Alphas[ 0 ] >= threshold ) )
      {
        // Remove the two points from the point set
        if ( meshCollection->GetPointParameters()->ElementAt( point0Id ).m_CanMoveX &&
             meshCollection->GetPointParameters()->ElementAt( point0Id ).m_CanMoveY &&
             meshCollection->GetPointParameters()->ElementAt( point0Id ).m_CanMoveZ )
        {
          sparsePoints->DeleteIndex( point0Id );
        }

        if ( meshCollection->GetPointParameters()->ElementAt( point1Id ).m_CanMoveX &&
             meshCollection->GetPointParameters()->ElementAt( point1Id ).m_CanMoveY &&
             meshCollection->GetPointParameters()->ElementAt( point1Id ).m_CanMoveZ )
        {
          sparsePoints->DeleteIndex( point1Id );
        }
      }

    }


    // Now write out a file that tetgen can read
    const std::string  tetgenInputFileName = "pointSet.node";
    std::ofstream  out( tetgenInputFileName.c_str() );
    if ( out.bad() )
    {
      std::cerr << "Can't open " << tetgenInputFileName << " for writing." << std::endl;
      return -1;
    }
    out << sparsePoints->Size() << " 3  0  0" << std::endl;
    int counter = 0;
    for ( kvl::AtlasMesh::PointsContainer::ConstIterator  it = sparsePoints->Begin();
          it != sparsePoints->End(); ++it, ++counter )
    {
      out << counter << "  " << it.Value()[ 0 ] << "  " << it.Value()[ 1 ] << "  " << it.Value()[ 2 ] << std::endl;
    }
    out.close();
    std::cout << "Wrote " << tetgenInputFileName  << std::endl;


    // Run tetgen
    std::cout << "Now do the following: \n" << std::endl;
    std::cout << "         ./tetgen " << tetgenInputFileName << std::endl;
    std::cout << "\nand then feed in any character to continue" << std::endl;
    char  dummy;
    std::cin >> dummy;


    // Re-index the sparse points to a dense representation
    kvl::AtlasMesh::PointsContainer::Pointer  reindexedSparsePoints = kvl::AtlasMesh::PointsContainer::New();
    kvl::AtlasMesh::PointIdentifier  newIndex = 0;
    for ( kvl::AtlasMesh::PointsContainer::ConstIterator  it = sparsePoints->Begin();
          it != sparsePoints->End(); ++it, ++newIndex )
    {
      reindexedSparsePoints->InsertElement( newIndex, it.Value() );
    }


    // Read the output from tetgen and make a mesh from it
    const std::string  tetgenOutputFileName = "pointSet.1.ele";
    std::ifstream  in( tetgenOutputFileName.c_str() );
    if ( in.bad() )
    {
      std::cerr << "Couldn't read from file " << tetgenOutputFileName << std::endl;
      return -1;
    }

    typedef itk::AutomaticTopologyMeshSource< kvl::AtlasMesh >  MeshSourceType;
    MeshSourceType::Pointer  meshSource = MeshSourceType::New();
    // const int size = 255;
    // char buffer[ size ];
    // const std::string  line = buffer;
    // std::istringstream  lineStream( line );
    // int  numberOfTetrahedra;
    // lineStream >> numberOfTetrahedra;
    // std::cout << "numberOfTetrahedra:" << numberOfTetrahedra << std::endl;
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

      // Get the corresponding points
      const kvl::AtlasMesh::PointType&  point0 = reindexedSparsePoints->ElementAt( point0Id );
      const kvl::AtlasMesh::PointType&  point1 = reindexedSparsePoints->ElementAt( point1Id );
      const kvl::AtlasMesh::PointType&  point2 = reindexedSparsePoints->ElementAt( point2Id );
      const kvl::AtlasMesh::PointType&  point3 = reindexedSparsePoints->ElementAt( point3Id );
      std::cout << "Adding tetrahedron with points " << std::endl;
      std::cout << "            " << point0 << std::endl;
      std::cout << "            " << point1 << std::endl;
      std::cout << "            " << point2 << std::endl;
      std::cout << "            " << point3 << std::endl;

      meshSource->AddTetrahedron( meshSource->AddPoint( point0 ),
                                  meshSource->AddPoint( point1 ),
                                  meshSource->AddPoint( point2 ),
                                  meshSource->AddPoint( point3 ) );
    }
    in.close();


    // Now complete the mesh produced by the mesh source with point parameters
    // Assign flat alphas as a starting point. Vertices lying on the border can not move freely and belong to
    // first class
    kvl::AtlasAlphasType   flatAlphasEntry( numberOfClasses );
    flatAlphasEntry.Fill( 1.0f / static_cast< float >( numberOfClasses ) );

    kvl::AtlasAlphasType   borderAlphasEntry( numberOfClasses );
    borderAlphasEntry.Fill( 0.0f );
    borderAlphasEntry[ 0 ] = 1.0f;

    for ( kvl::AtlasMesh::PointsContainer::ConstIterator  pointIt = meshSource->GetOutput()->GetPoints()->Begin();
          pointIt != meshSource->GetOutput()->GetPoints()->End();
          ++pointIt )
    {
      kvl::AtlasMesh::PixelType  pointParameters;

      pointParameters.m_Alphas = flatAlphasEntry;
      pointParameters.m_CanChangeAlphas = true;

      if ( ( pointIt.Value()[ 0 ] == 0 ) || ( pointIt.Value()[ 0 ] == ( domainSize[ 0 ] - 1 ) ) )
      {
        pointParameters.m_CanMoveX = false;
        pointParameters.m_Alphas = borderAlphasEntry;
        pointParameters.m_CanChangeAlphas = false;
      }
      else
      {
        pointParameters.m_CanMoveX = true;
      }

      if ( ( pointIt.Value()[ 1 ] == 0 ) || ( pointIt.Value()[ 1 ] == ( domainSize[ 1 ] - 1 ) ) )
      {
        pointParameters.m_CanMoveY = false;
        pointParameters.m_Alphas = borderAlphasEntry;
        pointParameters.m_CanChangeAlphas = false;
      }
      else
      {
        pointParameters.m_CanMoveY = true;
      }

      if ( ( pointIt.Value()[ 2 ] == 0 ) || ( pointIt.Value()[ 2 ] == ( domainSize[ 2 ] - 1 ) ) )
      {
        pointParameters.m_CanMoveZ = false;
        pointParameters.m_Alphas = borderAlphasEntry;
        pointParameters.m_CanChangeAlphas = false;
      }
      else
      {
        pointParameters.m_CanMoveZ = true;
      }

      meshSource->GetOutput()->SetPointData( pointIt.Index(), pointParameters );
    }


    // Now simply create the mesh collection by repeating the mesh
    const std::string  sparseMeshCollectionFileName = "sparseMeshCollection.txt";
    kvl::AtlasMeshCollection::Pointer  sparseMeshCollection = kvl::AtlasMeshCollection::New();
    sparseMeshCollection->GenerateFromSingleMesh( meshSource->GetOutput(), numberOfMeshes, stiffness );
    sparseMeshCollection->Write( sparseMeshCollectionFileName.c_str() );
    std::cout << "Wrote " << sparseMeshCollectionFileName << std::endl;


    // Estimate the model parameters for this mesh
    estimator->SetInitialMeshCollection( sparseMeshCollection );
    estimator->SetAlphasSmoothingFactor( 2 );
    estimator->Estimate();
    sparseMeshCollection->Write( sparseMeshCollectionFileName.c_str() );
    std::cout << "Wrote " << sparseMeshCollectionFileName << std::endl;

  }
  catch( itk::ExceptionObject& e )
  {
    std::cerr << e << std::endl;
  }

  return 0;
};

