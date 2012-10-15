/**
 * @file  kvlLabelImageToCompactMeshRepresentation.cxx
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
#include "kvlAtlasMeshCollection.h"
#include "kvlAtlasMeshCollectionValidator.h"
#include "kvlAtlasParameterEstimator.h"
#include "itkImage.h"
#include "itkImageRegionIteratorWithIndex.h"
#include "itkImageFileWriter.h"

#include "vnl/vnl_sample.h"


// Useful class for permutating edge list
class EdgeElement
{
public:
  //
  EdgeElement()
  {
    m_EdgeId = 0;
    m_EdgePriority = vnl_sample_uniform( 0,  1 );
  }

  //
  EdgeElement( kvl::AtlasMesh::CellIdentifier edgeId )
  {
    m_EdgeId = edgeId;
    m_EdgePriority = vnl_sample_uniform( 0,  1 );
  }

  //
  ~EdgeElement() {}

  // Needed for std::sort()
  bool operator<(const EdgeElement& it) const
  {
    return m_EdgePriority < it.m_EdgePriority;
  }

  //
  kvl::AtlasMesh::CellIdentifier  GetEdgeId() const
  {
    return m_EdgeId;
  }

private:
  kvl::AtlasMesh::CellIdentifier   m_EdgeId;
  double  m_EdgePriority;
};



std::vector< kvl::AtlasMesh::CellIdentifier >
Permutate(  const std::vector< kvl::AtlasMesh::CellIdentifier >&  edgeList )
{

  // Populate sortable vector
  std::vector< EdgeElement >    edgeElements;
  for ( std::vector< kvl::AtlasMesh::CellIdentifier >::const_iterator  it = edgeList.begin(); it != edgeList.end(); ++it )
  {
    edgeElements.push_back( EdgeElement( *it ) );
  }

  // Sort the sortable vector
  std::sort( edgeElements.begin(), edgeElements.end() );

  // Construct something to return
  std::vector< kvl::AtlasMesh::CellIdentifier >  result;
  for ( std::vector< EdgeElement >::const_iterator  it = edgeElements.begin(); it != edgeElements.end(); ++it )
  {
    result.push_back( ( *it ).GetEdgeId()  );
  }
  return result;
}



int main( int, char** )
{

  // Create an empty image
  typedef itk::Image< unsigned char, 3 >  ImageType;
  ImageType::SizeType  size;
  size[ 0 ] = 30;
  size[ 1 ] = 30;
  size[ 2 ] = 30;
  ImageType::Pointer  image = ImageType::New();
  image->SetRegions( size );
  image->Allocate();
  image->FillBuffer( 0 );

  // Draw a binary object in the image
  float  center[ 3 ];
  for ( int i = 0; i < 3; i++ )
  {
    center[ i ] = ( size[ i ] - 1 ) / 2.0f;
  }
  const float  radius = size[ 0 ] / 3.5f;
  for ( itk::ImageRegionIteratorWithIndex< ImageType >  it( image, image->GetBufferedRegion() );
        !it.IsAtEnd(); ++it )
  {
    float  distance = 0;
    for ( int i = 0; i < 3; i++ )
    {
      distance += pow( it.GetIndex()[ i ] - center[ i ], 2 );
    }
    distance = sqrt( distance );

    if ( distance < radius )
    {
      it.Value() = 1;
    }

  }

  typedef itk::ImageFileWriter< ImageType >  WriterType;
  WriterType::Pointer  writer = WriterType::New();
  writer->SetFileName( "inputImage.mhd" );
  writer->SetInput( image );
  writer->Update();


  // Construct a mesh collection, so that each node coincides exactly with each pixel
  // in the image
  kvl::AtlasMeshCollection::Pointer  collection =  kvl::AtlasMeshCollection::New();
  unsigned int  meshSize[ 3 ];
  unsigned int  domainSize[ 3 ];
  for ( int i = 0; i < 3; i++ )
  {
    meshSize[ i ] = size[ i ];
    domainSize[ i ] = size[ i ];
  }
  const unsigned int  numberOfClasses = 2;
  const unsigned int  numberOfMeshes = 1;
  collection->Construct( meshSize, domainSize, 1000,
                         numberOfClasses, numberOfMeshes );
  collection->Write( "initialMesh.txt" );


  // Estimate the model parameters for this mesh collection. This should yield probabilies
  // of either 0 or 1 in each of the nodes
  kvl::AtlasParameterEstimator::Pointer  estimator = kvl::AtlasParameterEstimator::New();
  std::vector< ImageType::ConstPointer >  labelImages;
  labelImages.push_back( image.GetPointer() );
  estimator->SetLabelImages( labelImages );
  estimator->SetInitialMeshCollection( collection );
  estimator->Estimate();
  //collection = estimator->GetCurrentMeshCollection();
  collection->Write( "estimatedMesh.txt" );



  // Now immobilize the points lying on the border. We need to do this so that their
  // position won't change during the edge collapse operation
  kvl::AtlasMesh::PointDataContainer::Pointer  originalPointParameters = kvl::AtlasMesh::PointDataContainer::New();
  for ( kvl::AtlasMesh::PointDataContainer::ConstIterator  pointParamIt = collection->GetPointParameters()->Begin();
        pointParamIt != collection->GetPointParameters()->End(); ++pointParamIt )
  {
    // Copy
    originalPointParameters->InsertElement( pointParamIt.Index(),  pointParamIt.Value() );
  }

  const kvl::AtlasMesh::CellLinksContainerPointer  cellLinks = collection->GetCellLinks();
  for ( kvl::AtlasMesh::PointDataContainer::Iterator  pointParamIt = collection->GetPointParameters()->Begin();
        pointParamIt != collection->GetPointParameters()->End(); ++pointParamIt )
  {

    // Loop over all the edges containing this point to check if it lies on a boundary
    bool  isBoundaryPoint = false;
    const std::set< kvl::AtlasMesh::CellIdentifier >&  cellsContainingPoint = cellLinks->ElementAt( pointParamIt.Index() );
    for ( std::set< kvl::AtlasMesh::CellIdentifier >::const_iterator  cellIt = cellsContainingPoint.begin();
          cellIt != cellsContainingPoint.end(); ++cellIt )
    {
      // Get the line cell
      const kvl::AtlasMesh::CellType*  cell = collection->GetCells()->ElementAt( *cellIt );
      if ( cell->GetType() != kvl::AtlasMesh::CellType::LINE_CELL )
      {
        continue;
      }

      // Retrieve ids of the points on the edge
      kvl::AtlasMesh::CellType::PointIdConstIterator  pointIt = cell->PointIdsBegin();
      kvl::AtlasMesh::PointIdentifier  point0Id = *pointIt;
      ++pointIt;
      kvl::AtlasMesh::PointIdentifier  point1Id = *pointIt;

      // Look up alphas in the point under investigation, and in the two edge points.
      // If there is a difference, we have an boundary point
      if ( ( pointParamIt.Value().m_Alphas != collection->GetPointParameters()->ElementAt( point0Id ).m_Alphas ) ||
           ( pointParamIt.Value().m_Alphas != collection->GetPointParameters()->ElementAt( point1Id ).m_Alphas ) )
      {
        isBoundaryPoint = true;
        break;
      }

    }

    //
    if ( isBoundaryPoint )
    {
      // Immobilize point
      pointParamIt.Value().m_CanMoveX = false;
      pointParamIt.Value().m_CanMoveY = false;
      pointParamIt.Value().m_CanMoveZ = false;
    }

  } // End loop over all points

  collection->Write( "immobilizedEstimatedMesh.txt" );



  // Walk several times over all the edges in the mesh, and try to collapse them if
  // appropriate and possible.
  for ( int iterationNumber = 0; ; iterationNumber++ )
  {
    int  numberOfEdgesCollapsed = 0;
    // Cache edges at this point
    std::vector< kvl::AtlasMesh::CellIdentifier >  edgesToTry;
    for ( kvl::AtlasMesh::CellsContainer::ConstIterator  cellIt = collection->GetCells()->Begin();
          cellIt != collection->GetCells()->End(); ++cellIt )
    {
      if ( cellIt.Value()->GetType() == kvl::AtlasMesh::CellType::LINE_CELL )
      {
        edgesToTry.push_back( cellIt.Index() );
      }
    }

    // Randomize the edges' order
    std::cout << "Randomizing edge order..." << std::endl;
    std::cout << "   number of edges before: " << edgesToTry.size() << std::endl;
    edgesToTry = Permutate(  edgesToTry );
    std::cout << "   number of edges after: " << edgesToTry.size() << std::endl;

    // Loop over all edges
    unsigned int  progressCounter = 0;
    float  progress = 0.0f;
    for ( std::vector< kvl::AtlasMesh::CellIdentifier >::const_iterator  edgeIt = edgesToTry.begin();
          edgeIt != edgesToTry.end(); ++edgeIt, ++progressCounter )
    {
      kvl::AtlasMesh::CellIdentifier  edgeId = *edgeIt;
      progress = static_cast< float >( progressCounter ) / static_cast< float >( edgesToTry.size() );
      std::cerr << "Progress: " << progress * 100 << "%" << std::endl;

      // Check if the edge is still present in the mesh.
      if ( !collection->GetCells()->IndexExists( edgeId ) )
      {
        continue;
      }


      std::cout << "    Analyzing edge with id: " << edgeId << std::endl;


      // Check if the two vertices of this edge have the same label
      const kvl::AtlasMesh::CellType*  edge = collection->GetCells()->ElementAt( edgeId );
      kvl::AtlasMesh::CellType::PointIdConstIterator  pointIt = edge->PointIdsBegin();
      kvl::AtlasMesh::PointIdentifier  edgePoint0Id = *pointIt;
      ++pointIt;
      kvl::AtlasMesh::PointIdentifier  edgePoint1Id = *pointIt;

      std::cout << " Comparing " << collection->GetPointParameters()->ElementAt( edgePoint0Id ).m_Alphas
                << " to " << collection->GetPointParameters()->ElementAt( edgePoint1Id ).m_Alphas << std::endl;
      if ( collection->GetPointParameters()->ElementAt( edgePoint0Id ).m_Alphas !=
           collection->GetPointParameters()->ElementAt( edgePoint1Id ).m_Alphas )
      {
        std::cout << "           Not the same; discarding this edge" << std::endl;
        continue;
      }
      else
      {
        std::cout  << "           The same!"<< std::endl;
      }


      // Because we're using such an ugly hack to do all this, points on the object boundary actually
      // would pull points on the mesh boundary inwards. Detect if that happens, and if it does, the
      // edge can't be collapsed
      /*      const  bool  point0LiesOnMeshBoundary =
                  ( ( collection->GetReferencePosition()->ElementAt( edgePoint0Id )[ 0 ] == 0 ) ||
                    ( collection->GetReferencePosition()->ElementAt( edgePoint0Id )[ 0 ] == ( size[ 0 ] - 1 ) ) ||
                    ( collection->GetReferencePosition()->ElementAt( edgePoint0Id )[ 1 ] == 0 ) ||
                    ( collection->GetReferencePosition()->ElementAt( edgePoint0Id )[ 1 ] == ( size[ 1 ] - 1 ) ) ||
                    ( collection->GetReferencePosition()->ElementAt( edgePoint0Id )[ 2 ] == 0 ) ||
                    ( collection->GetReferencePosition()->ElementAt( edgePoint0Id )[ 2 ] == ( size[ 2 ] - 1 ) ) );*/
      const  bool  point0LiesOnCorner =
        ( ( ( collection->GetReferencePosition()->ElementAt( edgePoint0Id )[ 0 ] == 0 ) &&
            ( collection->GetReferencePosition()->ElementAt( edgePoint0Id )[ 1 ] == 0 ) &&
            ( collection->GetReferencePosition()->ElementAt( edgePoint0Id )[ 2 ] == 0 ) ) ||
          ( ( collection->GetReferencePosition()->ElementAt( edgePoint0Id )[ 0 ] == size[ 0 ]-1 ) &&
            ( collection->GetReferencePosition()->ElementAt( edgePoint0Id )[ 1 ] == 0 ) &&
            ( collection->GetReferencePosition()->ElementAt( edgePoint0Id )[ 2 ] == 0 ) ) ||
          ( ( collection->GetReferencePosition()->ElementAt( edgePoint0Id )[ 0 ] == size[ 0 ]-1 ) &&
            ( collection->GetReferencePosition()->ElementAt( edgePoint0Id )[ 1 ] == size[ 1 ]-1 ) &&
            ( collection->GetReferencePosition()->ElementAt( edgePoint0Id )[ 2 ] == 0 ) ) ||
          ( ( collection->GetReferencePosition()->ElementAt( edgePoint0Id )[ 0 ] == 0 ) &&
            ( collection->GetReferencePosition()->ElementAt( edgePoint0Id )[ 1 ] == size[ 1 ]-1 ) &&
            ( collection->GetReferencePosition()->ElementAt( edgePoint0Id )[ 2 ] == 0 ) ) ||
          ( ( collection->GetReferencePosition()->ElementAt( edgePoint0Id )[ 0 ] == 0 ) &&
            ( collection->GetReferencePosition()->ElementAt( edgePoint0Id )[ 1 ] == 0 ) &&
            ( collection->GetReferencePosition()->ElementAt( edgePoint0Id )[ 2 ] == size[ 2 ]-1 ) ) ||
          ( ( collection->GetReferencePosition()->ElementAt( edgePoint0Id )[ 0 ] == size[ 0 ]-1 ) &&
            ( collection->GetReferencePosition()->ElementAt( edgePoint0Id )[ 1 ] == 0 ) &&
            ( collection->GetReferencePosition()->ElementAt( edgePoint0Id )[ 2 ] == size[ 2 ]-1 ) ) ||
          ( ( collection->GetReferencePosition()->ElementAt( edgePoint0Id )[ 0 ] == size[ 0 ]-1 ) &&
            ( collection->GetReferencePosition()->ElementAt( edgePoint0Id )[ 1 ] == size[ 1 ]-1 ) &&
            ( collection->GetReferencePosition()->ElementAt( edgePoint0Id )[ 2 ] == size[ 2 ]-1 ) ) ||
          ( ( collection->GetReferencePosition()->ElementAt( edgePoint0Id )[ 0 ] == 0 ) &&
            ( collection->GetReferencePosition()->ElementAt( edgePoint0Id )[ 1 ] == size[ 1 ]-1 ) &&
            ( collection->GetReferencePosition()->ElementAt( edgePoint0Id )[ 2 ] == size[ 2 ]-1 ) ) );

      int  numberOfConstraintsForPoint0 = 0;
      if ( !( collection->GetPointParameters()->ElementAt( edgePoint0Id ).m_CanMoveX ) )
      {
        numberOfConstraintsForPoint0++;
      }
      if ( !( collection->GetPointParameters()->ElementAt( edgePoint0Id ).m_CanMoveY ) )
      {
        numberOfConstraintsForPoint0++;
      }
      if ( !( collection->GetPointParameters()->ElementAt( edgePoint0Id ).m_CanMoveZ ) )
      {
        numberOfConstraintsForPoint0++;
      }
      const  bool  point0IsImmobilized = ( numberOfConstraintsForPoint0 == 3 );

      /*      const  bool  point1LiesOnMeshBoundary =
                  ( ( collection->GetReferencePosition()->ElementAt( edgePoint1Id )[ 0 ] == 0 ) ||
                    ( collection->GetReferencePosition()->ElementAt( edgePoint1Id )[ 0 ] == ( size[ 0 ] - 1 ) ) ||
                    ( collection->GetReferencePosition()->ElementAt( edgePoint1Id )[ 1 ] == 0 ) ||
                    ( collection->GetReferencePosition()->ElementAt( edgePoint1Id )[ 1 ] == ( size[ 1 ] - 1 ) ) ||
                    ( collection->GetReferencePosition()->ElementAt( edgePoint1Id )[ 2 ] == 0 ) ||
                    ( collection->GetReferencePosition()->ElementAt( edgePoint1Id )[ 2 ] == ( size[ 2 ] - 1 ) ) );*/
      const  bool  point1LiesOnCorner =
        ( ( ( collection->GetReferencePosition()->ElementAt( edgePoint1Id )[ 0 ] == 0 ) &&
            ( collection->GetReferencePosition()->ElementAt( edgePoint1Id )[ 1 ] == 0 ) &&
            ( collection->GetReferencePosition()->ElementAt( edgePoint1Id )[ 2 ] == 0 ) ) ||
          ( ( collection->GetReferencePosition()->ElementAt( edgePoint1Id )[ 0 ] == size[ 0 ]-1 ) &&
            ( collection->GetReferencePosition()->ElementAt( edgePoint1Id )[ 1 ] == 0 ) &&
            ( collection->GetReferencePosition()->ElementAt( edgePoint1Id )[ 2 ] == 0 ) ) ||
          ( ( collection->GetReferencePosition()->ElementAt( edgePoint1Id )[ 0 ] == size[ 0 ]-1 ) &&
            ( collection->GetReferencePosition()->ElementAt( edgePoint1Id )[ 1 ] == size[ 1 ]-1 ) &&
            ( collection->GetReferencePosition()->ElementAt( edgePoint1Id )[ 2 ] == 0 ) ) ||
          ( ( collection->GetReferencePosition()->ElementAt( edgePoint1Id )[ 0 ] == 0 ) &&
            ( collection->GetReferencePosition()->ElementAt( edgePoint1Id )[ 1 ] == size[ 1 ]-1 ) &&
            ( collection->GetReferencePosition()->ElementAt( edgePoint1Id )[ 2 ] == 0 ) ) ||
          ( ( collection->GetReferencePosition()->ElementAt( edgePoint1Id )[ 0 ] == 0 ) &&
            ( collection->GetReferencePosition()->ElementAt( edgePoint1Id )[ 1 ] == 0 ) &&
            ( collection->GetReferencePosition()->ElementAt( edgePoint1Id )[ 2 ] == size[ 2 ]-1 ) ) ||
          ( ( collection->GetReferencePosition()->ElementAt( edgePoint1Id )[ 0 ] == size[ 0 ]-1 ) &&
            ( collection->GetReferencePosition()->ElementAt( edgePoint1Id )[ 1 ] == 0 ) &&
            ( collection->GetReferencePosition()->ElementAt( edgePoint1Id )[ 2 ] == size[ 2 ]-1 ) ) ||
          ( ( collection->GetReferencePosition()->ElementAt( edgePoint1Id )[ 0 ] == size[ 0 ]-1 ) &&
            ( collection->GetReferencePosition()->ElementAt( edgePoint1Id )[ 1 ] == size[ 1 ]-1 ) &&
            ( collection->GetReferencePosition()->ElementAt( edgePoint1Id )[ 2 ] == size[ 2 ]-1 ) ) ||
          ( ( collection->GetReferencePosition()->ElementAt( edgePoint1Id )[ 0 ] == 0 ) &&
            ( collection->GetReferencePosition()->ElementAt( edgePoint1Id )[ 1 ] == size[ 1 ]-1 ) &&
            ( collection->GetReferencePosition()->ElementAt( edgePoint1Id )[ 2 ] == size[ 2 ]-1 ) ) );
      int  numberOfConstraintsForPoint1 = 0;
      if ( !( collection->GetPointParameters()->ElementAt( edgePoint1Id ).m_CanMoveX ) )
      {
        numberOfConstraintsForPoint1++;
      }
      if ( !( collection->GetPointParameters()->ElementAt( edgePoint1Id ).m_CanMoveY ) )
      {
        numberOfConstraintsForPoint1++;
      }
      if ( !( collection->GetPointParameters()->ElementAt( edgePoint1Id ).m_CanMoveZ ) )
      {
        numberOfConstraintsForPoint1++;
      }
      const  bool  point1IsImmobilized = ( numberOfConstraintsForPoint1 == 3 );

      const bool  immobilizedPoint0WantsToDrawInPoint1 =
        ( point0IsImmobilized && !point0LiesOnCorner &&
          ( numberOfConstraintsForPoint0 >=  numberOfConstraintsForPoint1 ) );

      const bool  immobilizedPoint1WantsToDrawInPoint0 =
        ( point1IsImmobilized && !point1LiesOnCorner &&
          ( numberOfConstraintsForPoint1 >=  numberOfConstraintsForPoint0 ) );

      std::cout << "Hack hack hack" << std::endl;
      std::cout << "       point0LiesOnCorner: " << point0LiesOnCorner << std::endl;
      std::cout << "       numberOfConstraintsForPoint0: " << numberOfConstraintsForPoint0 << std::endl;
      std::cout << "       point0IsImmobilized: " << point0IsImmobilized << std::endl;
      std::cout << std::endl;
      std::cout << "       point1LiesOnCorner: " << point1LiesOnCorner << std::endl;
      std::cout << "       numberOfConstraintsForPoint1: " << numberOfConstraintsForPoint1 << std::endl;
      std::cout << "       point1IsImmobilized: " << point1IsImmobilized << std::endl;
      std::cout << std::endl;
      std::cout << "       immobilizedPoint0WantsToDrawInPoint1: " <<  immobilizedPoint0WantsToDrawInPoint1 << std::endl;
      std::cout << "       immobilizedPoint1WantsToDrawInPoint0: " <<  immobilizedPoint1WantsToDrawInPoint0 << std::endl;

      if ( immobilizedPoint0WantsToDrawInPoint1 || immobilizedPoint1WantsToDrawInPoint0 )
      {
        // This is an artifically immobilized point. Check if the other point lies on mesh boundary
        std::cout << "           Artifically immobilized wants to pull mesh boundary point inwards. No way!" << std::endl;
        continue;
      }


      // Try to collapse
      kvl::AtlasMeshCollection::Pointer  collapsed;
      std::set< kvl::AtlasMesh::CellIdentifier >  disappearingCells;
      kvl::AtlasMesh::CellIdentifier  unifiedVertexId;
      if ( !collection->GetCollapsed( edgeId, collapsed, disappearingCells, unifiedVertexId, false ) )
      {
        continue;
      }



      // The resulting collapsed mesh is our new collection to work with
      numberOfEdgesCollapsed++;
      collection = collapsed;

      if ( size[ 0 ] < 6 )
      {
        // Write out
        std::ostringstream  fileNameStream;
        fileNameStream << "collection_" << iterationNumber << "_" << numberOfEdgesCollapsed
                       << "_collapsedEdgeId" << edgeId;
        collection->Write( fileNameStream.str().c_str() );

        // Validate collapsed edge mesh collection
        std::cout << "Validating collapsed" << std::endl;
        kvl::AtlasMeshCollectionValidator::Pointer  validator = kvl::AtlasMeshCollectionValidator::New();
        if ( !validator->Validate( collapsed ) )
        {
          exit( -1 );
        }
      }

      //char  input;
      //std::cout << "Type character to continue" << std::endl;
      //std::cin >> input;
    } // End loop over all edges to be analyzed


    // Write out what we have now
    std::ostringstream  fileNameStream;
    fileNameStream << "collection_" << iterationNumber;
    collection->Write( fileNameStream.str().c_str() );


    // Check if we can stop
    if ( numberOfEdgesCollapsed == 0 )
    {
      std::cout << "No more edges to try; stopping. :-)" << std::endl;
      break;
    }
  } // End repeated interations


  // Restore original mobility of the points
  for ( kvl::AtlasMesh::PointDataContainer::Iterator  pointParamIt = collection->GetPointParameters()->Begin();
        pointParamIt != collection->GetPointParameters()->End(); ++pointParamIt )
  {
    // Copy
    pointParamIt.Value().m_CanMoveX = originalPointParameters->ElementAt( pointParamIt.Index() ).m_CanMoveX;
    pointParamIt.Value().m_CanMoveY = originalPointParameters->ElementAt( pointParamIt.Index() ).m_CanMoveY;
    pointParamIt.Value().m_CanMoveZ = originalPointParameters->ElementAt( pointParamIt.Index() ).m_CanMoveZ;
  }



  return 0;
};

