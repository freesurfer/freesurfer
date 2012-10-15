/**
 * @file  kvlAtlasMeshVertexProcessor.cxx
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
#include "kvlAtlasMeshVertexProcessor.h"

#include "kvlAtlasMeshRasterizor.h"

#if 1
#include "kvlAtlasMeshPositionGradientCalculator.h"
#endif

namespace kvl
{


//
//
//
AtlasMeshVertexProcessor
::AtlasMeshVertexProcessor()
{
  m_MeshCollection = 0;
  m_Beta = 0.0f;
}



//
//
//
AtlasMeshVertexProcessor
::~AtlasMeshVertexProcessor()
{
}




//
//
//
void
AtlasMeshVertexProcessor
::PrintSelf( std::ostream& os, itk::Indent indent ) const
{
}


//
//
//
void
AtlasMeshVertexProcessor
::SetLabelImages( const std::vector< LabelImageType::ConstPointer >& labelImages )
{
  m_LabelImages = labelImages;
}


//
//
//
void
AtlasMeshVertexProcessor
::
CalculateVertexNeighborhoods()
{
  // Sanity check
  if ( !m_MeshCollection )
  {
    itkExceptionMacro( "No mesh collection!" );
  }


  // Initialize each vertex's neighboring triangles to an empty container
  m_VertexNeighborhoodsContainers.clear();
  for ( int meshNumber = 0; meshNumber < m_MeshCollection->GetNumberOfMeshes(); meshNumber++ )
  {
    VertexNeighborhoodsContainerType::Pointer  vertexNeighborhoods = VertexNeighborhoodsContainerType::New();

    AtlasMesh::PointsContainer::ConstIterator  pointIt = m_MeshCollection->GetPositions()[ 0 ]->Begin();
    while ( pointIt != m_MeshCollection->GetPositions()[ 0 ]->End() )
    {
      vertexNeighborhoods->CreateIndex( pointIt.Index() );

      ++pointIt;
    }

    m_VertexNeighborhoodsContainers.push_back( vertexNeighborhoods );
  }


  // Loop over all tetrahedra
  for ( AtlasMesh::CellsContainer::ConstIterator  cellIt = m_MeshCollection->GetCells()->Begin();
        cellIt != m_MeshCollection->GetCells()->End(); ++cellIt )
  {
    const AtlasMesh::CellType*  cell = cellIt.Value();

    if( cell->GetType() != AtlasMesh::CellType::TETRAHEDRON_CELL )
    {
      continue;
    }

    // Retrieve point id's of the four vertices
    AtlasMesh::CellType::PointIdConstIterator pointIt = cell->PointIdsBegin();
    const AtlasMesh::PointIdentifier point0Id = *pointIt;
    ++pointIt;
    const AtlasMesh::PointIdentifier point1Id = *pointIt;
    ++pointIt;
    const AtlasMesh::PointIdentifier point2Id = *pointIt;
    ++pointIt;
    const AtlasMesh::PointIdentifier point3Id = *pointIt;

    // Retrieve address of tetrahedron info stuff
    const float*  referenceVolumeTimesK =
      &( m_MeshCollection->GetReferenceTetrahedronInfos()->ElementAt( cellIt.Index() ).m_ReferenceVolumeTimesK );

    const float*  z11 =
      &( m_MeshCollection->GetReferenceTetrahedronInfos()->ElementAt( cellIt.Index() ).m_Z11 );
    const float*  z21 =
      &( m_MeshCollection->GetReferenceTetrahedronInfos()->ElementAt( cellIt.Index() ).m_Z21 );
    const float*  z31 =
      &( m_MeshCollection->GetReferenceTetrahedronInfos()->ElementAt( cellIt.Index() ).m_Z31 );
    const float*  z41 =
      &( m_MeshCollection->GetReferenceTetrahedronInfos()->ElementAt( cellIt.Index() ).m_Z41 );

    const float*  z12 =
      &( m_MeshCollection->GetReferenceTetrahedronInfos()->ElementAt( cellIt.Index() ).m_Z12 );
    const float*  z22 =
      &( m_MeshCollection->GetReferenceTetrahedronInfos()->ElementAt( cellIt.Index() ).m_Z22 );
    const float*  z32 =
      &( m_MeshCollection->GetReferenceTetrahedronInfos()->ElementAt( cellIt.Index() ).m_Z32 );
    const float*  z42 =
      &( m_MeshCollection->GetReferenceTetrahedronInfos()->ElementAt( cellIt.Index() ).m_Z42 );

    const float*  z13 =
      &( m_MeshCollection->GetReferenceTetrahedronInfos()->ElementAt( cellIt.Index() ).m_Z13 );
    const float*  z23 =
      &( m_MeshCollection->GetReferenceTetrahedronInfos()->ElementAt( cellIt.Index() ).m_Z23 );
    const float*  z33 =
      &( m_MeshCollection->GetReferenceTetrahedronInfos()->ElementAt( cellIt.Index() ).m_Z33 );
    const float*  z43 =
      &( m_MeshCollection->GetReferenceTetrahedronInfos()->ElementAt( cellIt.Index() ).m_Z43 );


    // Loop over all meshes
    for ( int meshNumber = 0; meshNumber < m_MeshCollection->GetNumberOfMeshes(); meshNumber++ )
    {
      // Retrieve address of (x, y, z) values in each of the vertices
      AtlasMesh::PointsContainer::Pointer positions = m_MeshCollection->GetPositions()[ meshNumber ];
      const AtlasMesh::PointType::ValueType*  x0 = &( positions->ElementAt( point0Id ).GetElement( 0 ) );
      const AtlasMesh::PointType::ValueType*  y0 = &( positions->ElementAt( point0Id ).GetElement( 1 ) );
      const AtlasMesh::PointType::ValueType*  z0 = &( positions->ElementAt( point0Id ).GetElement( 2 ) );
      const AtlasMesh::PointType::ValueType*  x1 = &( positions->ElementAt( point1Id ).GetElement( 0 ) );
      const AtlasMesh::PointType::ValueType*  y1 = &( positions->ElementAt( point1Id ).GetElement( 1 ) );
      const AtlasMesh::PointType::ValueType*  z1 = &( positions->ElementAt( point1Id ).GetElement( 2 ) );
      const AtlasMesh::PointType::ValueType*  x2 = &( positions->ElementAt( point2Id ).GetElement( 0 ) );
      const AtlasMesh::PointType::ValueType*  y2 = &( positions->ElementAt( point2Id ).GetElement( 1 ) );
      const AtlasMesh::PointType::ValueType*  z2 = &( positions->ElementAt( point2Id ).GetElement( 2 ) );
      const AtlasMesh::PointType::ValueType*  x3 = &( positions->ElementAt( point3Id ).GetElement( 0 ) );
      const AtlasMesh::PointType::ValueType*  y3 = &( positions->ElementAt( point3Id ).GetElement( 1 ) );
      const AtlasMesh::PointType::ValueType*  z3 = &( positions->ElementAt( point3Id ).GetElement( 2 ) );

      // Retrieve address of alphas in each of the vertices
      const AtlasAlphasType*  alphas0 =
        &( m_MeshCollection->GetPointParameters()->ElementAt( point0Id ).m_Alphas );
      const AtlasAlphasType*  alphas1 =
        &( m_MeshCollection->GetPointParameters()->ElementAt( point1Id ).m_Alphas );
      const AtlasAlphasType*  alphas2 =
        &( m_MeshCollection->GetPointParameters()->ElementAt( point2Id ).m_Alphas );
      const AtlasAlphasType*  alphas3 =
        &( m_MeshCollection->GetPointParameters()->ElementAt( point3Id ).m_Alphas );

      // Add this tetrahedron to the neighborhood of vertex 0
      VertexNeighboringTetrahedronInfo  vertex0Info;
      vertex0Info.m_TetrahedronId = cellIt.Index();
      vertex0Info.m_ReferenceVolumeTimesK = referenceVolumeTimesK;
      vertex0Info.m_Z11 = z11;
      vertex0Info.m_Z21 = z21;
      vertex0Info.m_Z31 = z31;
      vertex0Info.m_Z41 = z41;
      vertex0Info.m_Z12 = z12;
      vertex0Info.m_Z22 = z22;
      vertex0Info.m_Z32 = z32;
      vertex0Info.m_Z42 = z42;
      vertex0Info.m_Z13 = z13;
      vertex0Info.m_Z23 = z23;
      vertex0Info.m_Z33 = z33;
      vertex0Info.m_Z43 = z43;
      vertex0Info.m_X1 = x1;
      vertex0Info.m_Y1 = y1;
      vertex0Info.m_Z1 = z1;
      vertex0Info.m_X2 = x2;
      vertex0Info.m_Y2 = y2;
      vertex0Info.m_Z2 = z2;
      vertex0Info.m_X3 = x3;
      vertex0Info.m_Y3 = y3;
      vertex0Info.m_Z3 = z3;
      vertex0Info.m_Alphas0 = alphas0;
      vertex0Info.m_Alphas1 = alphas1;
      vertex0Info.m_Alphas2 = alphas2;
      vertex0Info.m_Alphas3 = alphas3;
      m_VertexNeighborhoodsContainers[ meshNumber ]->ElementAt( point0Id ).push_back( vertex0Info );

#if 0
      continue;
#endif

      // Add this tetrahedron to the neighborhood of vertex 1. We renumber the vertices based on the
      // following convention (just making sure we keep a "positive volume" tetrahedron):
      //    * newP0 = p1
      //    * newP1 = p2
      //    * newP2 = p0
      //    * newP3 = p3
      // The inverse of the matrix [newP0 newP1 newP2 newP3; 1 1 1 1 ] (elements of the "newZ") are
      // permutations of the original matrix Z. The exact permutation can be calculated by writing
      // newP = P * permutationMatrixWithZerosAndOnes, and then using
      // the fact that newZ = inv( newP )
      //                    = inv( permutationMatrixWithZerosAndOnes ) * inv( P )
      //                    = inv( permutationMatrixWithZerosAndOnes ) * Z
      //                    = otherPermutationMatrixWithZerosAndOnes * Z
      VertexNeighboringTetrahedronInfo  vertex1Info;
      vertex1Info.m_TetrahedronId = cellIt.Index();
      vertex1Info.m_ReferenceVolumeTimesK = referenceVolumeTimesK;
      vertex1Info.m_Z11 = z21;
      vertex1Info.m_Z21 = z31;
      vertex1Info.m_Z31 = z11;
      vertex1Info.m_Z41 = z41;
      vertex1Info.m_Z12 = z22;
      vertex1Info.m_Z22 = z32;
      vertex1Info.m_Z32 = z12;
      vertex1Info.m_Z42 = z42;
      vertex1Info.m_Z13 = z23;
      vertex1Info.m_Z23 = z33;
      vertex1Info.m_Z33 = z13;
      vertex1Info.m_Z43 = z43;
      vertex1Info.m_X1 = x2;
      vertex1Info.m_Y1 = y2;
      vertex1Info.m_Z1 = z2;
      vertex1Info.m_X2 = x0;
      vertex1Info.m_Y2 = y0;
      vertex1Info.m_Z2 = z0;
      vertex1Info.m_X3 = x3;
      vertex1Info.m_Y3 = y3;
      vertex1Info.m_Z3 = z3;
      vertex1Info.m_Alphas0 = alphas1;
      vertex1Info.m_Alphas1 = alphas2;
      vertex1Info.m_Alphas2 = alphas0;
      vertex1Info.m_Alphas3 = alphas3;
      m_VertexNeighborhoodsContainers[ meshNumber ]->ElementAt( point1Id ).push_back( vertex1Info );


      // Add this tetrahedron to the neighborhood of vertex 2. We renumber the vertices based on the
      // following convention (just making sure we keep a "positive volume" tetrahedron):
      //    * newP0 = p2
      //    * newP1 = p0
      //    * newP2 = p1
      //    * newP3 = p3
      // The inverse of the matrix [newP0 newP1 newP2 newP3; 1 1 1 1 ] (elements of the "newZ") are
      // permutations of the original matrix Z. The exact permutation can be calculated by writing
      // newP = P * permutationMatrixWithZerosAndOnes, and then using
      // the fact that newZ = inv( newP )
      //                    = inv( permutationMatrixWithZerosAndOnes ) * inv( P )
      //                    = inv( permutationMatrixWithZerosAndOnes ) * Z
      //                    = otherPermutationMatrixWithZerosAndOnes * Z
      VertexNeighboringTetrahedronInfo  vertex2Info;
      vertex2Info.m_TetrahedronId = cellIt.Index();
      vertex2Info.m_ReferenceVolumeTimesK = referenceVolumeTimesK;
      vertex2Info.m_Z11 = z31;
      vertex2Info.m_Z21 = z11;
      vertex2Info.m_Z31 = z21;
      vertex2Info.m_Z41 = z41;
      vertex2Info.m_Z12 = z32;
      vertex2Info.m_Z22 = z12;
      vertex2Info.m_Z32 = z22;
      vertex2Info.m_Z42 = z42;
      vertex2Info.m_Z13 = z33;
      vertex2Info.m_Z23 = z13;
      vertex2Info.m_Z33 = z23;
      vertex2Info.m_Z43 = z43;
      vertex2Info.m_X1 = x0;
      vertex2Info.m_Y1 = y0;
      vertex2Info.m_Z1 = z0;
      vertex2Info.m_X2 = x1;
      vertex2Info.m_Y2 = y1;
      vertex2Info.m_Z2 = z1;
      vertex2Info.m_X3 = x3;
      vertex2Info.m_Y3 = y3;
      vertex2Info.m_Z3 = z3;
      vertex2Info.m_Alphas0 = alphas2;
      vertex2Info.m_Alphas1 = alphas0;
      vertex2Info.m_Alphas2 = alphas1;
      vertex2Info.m_Alphas3 = alphas3;
      m_VertexNeighborhoodsContainers[ meshNumber ]->ElementAt( point2Id ).push_back( vertex2Info );


      // Add this tetrahedron to the neighborhood of vertex 3. We renumber the vertices based on the
      // following convention (just making sure we keep a "positive volume" tetrahedron):
      //    * newP0 = p3
      //    * newP1 = p1
      //    * newP2 = p0
      //    * newP3 = p2
      // The inverse of the matrix [newP0 newP1 newP2 newP3; 1 1 1 1 ] (elements of the "newZ") are
      // permutations of the original matrix Z. The exact permutation can be calculated by writing
      // newP = P * permutationMatrixWithZerosAndOnes, and then using
      // the fact that newZ = inv( newP )
      //                    = inv( permutationMatrixWithZerosAndOnes ) * inv( P )
      //                    = inv( permutationMatrixWithZerosAndOnes ) * Z
      //                    = otherPermutationMatrixWithZerosAndOnes * Z
      VertexNeighboringTetrahedronInfo  vertex3Info;
      vertex3Info.m_TetrahedronId = cellIt.Index();
      vertex3Info.m_ReferenceVolumeTimesK = referenceVolumeTimesK;
      vertex3Info.m_Z11 = z41;
      vertex3Info.m_Z21 = z21;
      vertex3Info.m_Z31 = z11;
      vertex3Info.m_Z41 = z31;
      vertex3Info.m_Z12 = z42;
      vertex3Info.m_Z22 = z22;
      vertex3Info.m_Z32 = z12;
      vertex3Info.m_Z42 = z32;
      vertex3Info.m_Z13 = z43;
      vertex3Info.m_Z23 = z23;
      vertex3Info.m_Z33 = z13;
      vertex3Info.m_Z43 = z33;
      vertex3Info.m_X1 = x1;
      vertex3Info.m_Y1 = y1;
      vertex3Info.m_Z1 = z1;
      vertex3Info.m_X2 = x0;
      vertex3Info.m_Y2 = y0;
      vertex3Info.m_Z2 = z0;
      vertex3Info.m_X3 = x2;
      vertex3Info.m_Y3 = y2;
      vertex3Info.m_Z3 = z2;
      vertex3Info.m_Alphas0 = alphas3;
      vertex3Info.m_Alphas1 = alphas1;
      vertex3Info.m_Alphas2 = alphas0;
      vertex3Info.m_Alphas3 = alphas2;
      m_VertexNeighborhoodsContainers[ meshNumber ]->ElementAt( point3Id ).push_back( vertex3Info );



    } // End loop over all meshes


  } // End loop over all tetrahedra



}



//
//
//
float
AtlasMeshVertexProcessor
::CalculateCost( int meshNumber, AtlasMesh::PointIdentifier pointId, float x, float y, float z )
{

  // Calculate prior contribution
  float  priorCost = this->CalculatePriorCost( meshNumber, pointId, x, y, z );

  if ( ( m_Beta == 0 ) || ( priorCost == itk::NumericTraits< float >::max() ) )
  {
    return priorCost;
  }

  if ( meshNumber >= m_LabelImages.size() )
  {
    itkExceptionMacro( << "Cannot determine data contribution to cost because the necessary label image"
                       " is not set" );
  }

  // Retrieve the neighborhood of this vertex
  const VertexNeighborhood&   neighborhood = m_VertexNeighborhoodsContainers[ meshNumber ]->ElementAt( pointId );

  // Calculate data contribution by looping over all tetrahedra and rendering each one
  typedef AtlasMeshRasterizor< FragmentProcessor::CalculateDataContributions >  DataContributionsCalculator;
  DataContributionsCalculator::Pointer  calculator = DataContributionsCalculator::New();
  calculator->SetLabelImage( m_LabelImages[ meshNumber ] );
  for ( VertexNeighborhood::const_iterator tetrahedronInfoIt = neighborhood.begin();
        tetrahedronInfoIt != neighborhood.end();
        ++tetrahedronInfoIt )
  {
    // Retrieve cached information about this tetrahedron
    const float x1 = *( ( *tetrahedronInfoIt ).m_X1 );
    const float y1 = *( ( *tetrahedronInfoIt ).m_Y1 );
    const float z1 = *( ( *tetrahedronInfoIt ).m_Z1 );
    const float x2 = *( ( *tetrahedronInfoIt ).m_X2 );
    const float y2 = *( ( *tetrahedronInfoIt ).m_Y2 );
    const float z2 = *( ( *tetrahedronInfoIt ).m_Z2 );
    const float x3 = *( ( *tetrahedronInfoIt ).m_X3 );
    const float y3 = *( ( *tetrahedronInfoIt ).m_Y3 );
    const float z3 = *( ( *tetrahedronInfoIt ).m_Z3 );
    const AtlasAlphasType*  alphas0 = ( *tetrahedronInfoIt ).m_Alphas0;
    const AtlasAlphasType*  alphas1 = ( *tetrahedronInfoIt ).m_Alphas1;
    const AtlasAlphasType*  alphas2 = ( *tetrahedronInfoIt ).m_Alphas2;
    const AtlasAlphasType*  alphas3 = ( *tetrahedronInfoIt ).m_Alphas3;

    // Rasterize the tetrahedron
    calculator->GetFragmentProcessor().StartNewTetrahedron( x, y, z,
        x1, y1, z1,
        x2, y2, z2,
        x3, y3, z3,
        *alphas0, *alphas1, *alphas2, *alphas3 );
    const float  pArray[] = {x, y, z };
    const AtlasMesh::PointType  p = pArray;
    const float  p1Array[] = {x1, y1, z1 };
    const AtlasMesh::PointType  p1 = p1Array;
    const float  p2Array[] = {x2, y2, z2 };
    const AtlasMesh::PointType  p2 = p2Array;
    const float  p3Array[] = {x3, y3, z3 };
    const AtlasMesh::PointType  p3 = p3Array;
    calculator->RasterizeTetrahedron( p, p1, p2, p3 );

    //std::cout << "Cost after rasterizing tetrahedron with id " << tetrahedronInfoIt->m_TetrahedronId
    //          << ": " << calculator->GetFragmentProcessor().GetCost() << std::endl;

  } // End loop over tetrahedra to calculate cost

  const float  dataCost = m_Beta * calculator->GetFragmentProcessor().GetCost();
  //std::cout << "dataCost: " << dataCost << std::endl;
  //std::cout << "priorCost: " << priorCost << std::endl;

  // Return sum of two contributions
  return dataCost + priorCost;
}




//
//
//
AtlasPositionGradientType
AtlasMeshVertexProcessor
::CalculateGradient( int meshNumber, AtlasMesh::PointIdentifier pointId, float x, float y, float z )
{
  // Calculate prior contribution
  AtlasPositionGradientType  priorGradient = this->CalculatePriorGradient( meshNumber, pointId, x, y, z );

  if ( m_Beta == 0 )
  {
    return priorGradient;
  }

  // Sanity checking
  if ( meshNumber >= m_LabelImages.size() )
  {
    itkExceptionMacro( << "Cannot determine data contribution to cost because the necessary label image"
                       " is not set" );
  }

  // Retrieve the neighborhood of this vertex
  const VertexNeighborhood&   neighborhood = m_VertexNeighborhoodsContainers[ meshNumber ]->ElementAt( pointId );

  // Calculate data contribution by looping over all tetrahedra and rendering each one
  typedef AtlasMeshRasterizor< FragmentProcessor::CalculateDataContributions >  DataContributionsCalculator;
  DataContributionsCalculator::Pointer  calculator = DataContributionsCalculator::New();
  calculator->SetLabelImage( m_LabelImages[ meshNumber ] );
  for ( VertexNeighborhood::const_iterator tetrahedronInfoIt = neighborhood.begin();
        tetrahedronInfoIt != neighborhood.end();
        ++tetrahedronInfoIt )
  {
    // Retrieve cached information about this tetrahedron
    const float x1 = *( ( *tetrahedronInfoIt ).m_X1 );
    const float y1 = *( ( *tetrahedronInfoIt ).m_Y1 );
    const float z1 = *( ( *tetrahedronInfoIt ).m_Z1 );
    const float x2 = *( ( *tetrahedronInfoIt ).m_X2 );
    const float y2 = *( ( *tetrahedronInfoIt ).m_Y2 );
    const float z2 = *( ( *tetrahedronInfoIt ).m_Z2 );
    const float x3 = *( ( *tetrahedronInfoIt ).m_X3 );
    const float y3 = *( ( *tetrahedronInfoIt ).m_Y3 );
    const float z3 = *( ( *tetrahedronInfoIt ).m_Z3 );
    const AtlasAlphasType*  alphas0 = ( *tetrahedronInfoIt ).m_Alphas0;
    const AtlasAlphasType*  alphas1 = ( *tetrahedronInfoIt ).m_Alphas1;
    const AtlasAlphasType*  alphas2 = ( *tetrahedronInfoIt ).m_Alphas2;
    const AtlasAlphasType*  alphas3 = ( *tetrahedronInfoIt ).m_Alphas3;

    // Rasterize the tetrahedron
    calculator->GetFragmentProcessor().StartNewTetrahedron( x, y, z,
        x1, y1, z1,
        x2, y2, z2,
        x3, y3, z3,
        *alphas0, *alphas1, *alphas2, *alphas3 );
    const float  pArray[] = {x, y, z };
    const AtlasMesh::PointType  p = pArray;
    const float  p1Array[] = {x1, y1, z1 };
    const AtlasMesh::PointType  p1 = p1Array;
    const float  p2Array[] = {x2, y2, z2 };
    const AtlasMesh::PointType  p2 = p2Array;
    const float  p3Array[] = {x3, y3, z3 };
    const AtlasMesh::PointType  p3 = p3Array;
    calculator->RasterizeTetrahedron( p, p1, p2, p3 );

  } // End loop over tetrahedra to calculate cost

  AtlasPositionGradientType  dataGradient = calculator->GetFragmentProcessor().GetGradient();
  dataGradient[ 0 ] *= m_Beta;
  dataGradient[ 1 ] *= m_Beta;
  dataGradient[ 2 ] *= m_Beta;

#if 1
  {
    std::cout << "CalculateGradient: dataGradient " << dataGradient << std::endl;
    std::cout << "CalculateGradient: piorGradient " << priorGradient << std::endl;
  }
#endif

  // Return sum of two contributions
  AtlasPositionGradientType  gradient;
  for ( int i = 0; i < 3; i++ )
  {
    gradient[ i ] = priorGradient[ i ];
    if ( gradient[ i ] != 0 )
    {
      gradient[ i ] += dataGradient[ i ];
    }
  }

  return gradient;
}


#if 0

//
//
//
Curvature
AtlasMeshVertexProcessor
::CalculateCurvature( int meshNumber, AtlasMesh::PointIdentifier pointId, float x, float y )
{
  // Calculate prior contribution
  Curvature  priorCurvature = this->CalculatePriorCurvature( meshNumber, pointId, x, y );

  if ( m_Beta == 0 )
  {
    return priorCurvature;
  }

  //
  if ( meshNumber >= m_LabelImages.size() )
  {
    itkExceptionMacro( << "Cannot determine data contribution to cost because the necessary label image"
                       " is not set" );
  }

  // Retrieve the neighborhood of this vertex
  const VertexNeighborhood&   neighborhood = m_VertexNeighborhoodsContainers[ meshNumber ]->ElementAt( pointId );

  // Calculate data contribution by looping over all triangles and rendering each one
  typedef AtlasMeshRasterizor< FragmentProcessor::CalculateDataContributions >  DataContributionsCalculator;
  DataContributionsCalculator::Pointer  calculator = DataContributionsCalculator::New();
  calculator->SetLabelImage( m_LabelImages[ meshNumber ] );
  for ( VertexNeighborhood::const_iterator triangleInfoIt = neighborhood.begin();
        triangleInfoIt != neighborhood.end();
        ++triangleInfoIt )
  {
    // Retrieve cached information about this triangle
    const float x1 = *( ( *triangleInfoIt ).m_X1 );
    const float y1 = *( ( *triangleInfoIt ).m_Y1 );
    const float x2 = *( ( *triangleInfoIt ).m_X2 );
    const float y2 = *( ( *triangleInfoIt ).m_Y2 );
    const AtlasAlphasType*  alphas0 = ( *triangleInfoIt ).m_Alphas0;
    const AtlasAlphasType*  alphas1 = ( *triangleInfoIt ).m_Alphas1;
    const AtlasAlphasType*  alphas2 = ( *triangleInfoIt ).m_Alphas2;

    // Rasterize the triangle
    calculator->GetFragmentProcessor().StartNewTriangle( x, y, x1, y1, x2, y2, *alphas0, *alphas1, *alphas2 );
    const float  p0[2] = { x, y };
    const float  p1[2] = { x1, y1 };
    const float  p2[2] = { x2, y2 };
    calculator->Rasterize( p0, p1, p2 );
  } // End loop over triangles to calculate cost

  Curvature  dataCurvature = calculator->GetFragmentProcessor().GetCurvature();
  dataCurvature.m_Curvature_dxdx *= m_Beta;
  dataCurvature.m_Curvature_dxdy *= m_Beta;
  dataCurvature.m_Curvature_dydy *= m_Beta;

  // Return sum of data and prior contributions
  Curvature  curvature;
  curvature.m_Curvature_dxdx = priorCurvature.m_Curvature_dxdx;
  curvature.m_Curvature_dxdy = priorCurvature.m_Curvature_dxdy;
  curvature.m_Curvature_dydy = priorCurvature.m_Curvature_dydy;
  if ( curvature.m_Curvature_dxdx != itk::NumericTraits< float >::max() )
  {
    curvature.m_Curvature_dxdx += dataCurvature.m_Curvature_dxdx;
    if ( curvature.m_Curvature_dydy != itk::NumericTraits< float >::max() )
    {
      curvature.m_Curvature_dxdy += dataCurvature.m_Curvature_dxdy;
    }
  }
  if ( curvature.m_Curvature_dydy != itk::NumericTraits< float >::max() )
  {
    curvature.m_Curvature_dydy += dataCurvature.m_Curvature_dydy;
  }

  return curvature;
}


#else

//
//
//
Curvature
AtlasMeshVertexProcessor
::CalculateCurvature( int meshNumber, AtlasMesh::PointIdentifier pointId, float x, float y, float z, bool verbose )
{
  // Sanity check on input
  if ( meshNumber >= m_LabelImages.size() )
  {
    itkExceptionMacro( << "Cannot determine data contribution to cost because the necessary label image"
                       " is not set" );
  }


  // Quick and dirty: use finite differences
  const float  delta = /* 0.005 */ /* 0.1 */ 2.0;
  float  xdelta = delta;
  float  ydelta = delta;
  float  zdelta = delta;
  if ( !m_MeshCollection->GetPointParameters()->ElementAt( pointId ).m_CanMoveX )
  {
    xdelta = 0;
  }
  if ( !m_MeshCollection->GetPointParameters()->ElementAt( pointId ).m_CanMoveY )
  {
    ydelta = 0;
  }
  if ( !m_MeshCollection->GetPointParameters()->ElementAt( pointId ).m_CanMoveZ )
  {
    zdelta = 0;
  }



  const float  cost = this->CalculateCost( meshNumber, pointId, x, y, z );
  const float  costPlusX = this->CalculateCost( meshNumber, pointId, x+xdelta, y, z );
  const float  costMinX = this->CalculateCost( meshNumber, pointId, x-xdelta, y, z );
  const float  costPlusY = this->CalculateCost( meshNumber, pointId, x, y+ydelta, z );
  const float  costMinY = this->CalculateCost( meshNumber, pointId, x, y-ydelta, z );
  const float  costPlusZ = this->CalculateCost( meshNumber, pointId, x, y, z+zdelta );
  const float  costMinZ = this->CalculateCost( meshNumber, pointId, x, y, z-zdelta );
  const float  costPlusXPlusY = this->CalculateCost( meshNumber, pointId, x+xdelta, y+ydelta, z );
  const float  costPlusXPlusZ = this->CalculateCost( meshNumber, pointId, x+xdelta, y, z+zdelta );
  const float  costPlusYPlusZ = this->CalculateCost( meshNumber, pointId, x, y+ydelta, z+zdelta );

  if ( verbose )
  {
    std::cout << "cost: " << cost << std::endl;
    std::cout << "costPlusX - cost: " << costPlusX - cost << std::endl;
    std::cout << "costMinX - cost: " << costMinX - cost << std::endl;
    std::cout << "costPlusY - cost: " << costPlusY - cost << std::endl;
    std::cout << "costMinY - cost: " << costMinY - cost << std::endl;
    std::cout << "costPlusZ - cost: " << costPlusZ - cost << std::endl;
    std::cout << "costMinZ - cost: " << costMinZ - cost << std::endl;
    std::cout << "costPlusXPlusY - cost: " << costPlusXPlusY - cost << std::endl;
    std::cout << "costPlusXPlusZ - cost: " << costPlusXPlusZ - cost << std::endl;
    std::cout << "costPlusYPlusZ - cost: " << costPlusYPlusZ - cost << std::endl;
  }


  Curvature  curvature;
  curvature.m_Curvature_dxdx = ( costPlusX - 2*cost + costMinX ) / ( delta * delta );
  curvature.m_Curvature_dydy = ( costPlusY - 2*cost + costMinY ) / ( delta * delta );
  curvature.m_Curvature_dzdz = ( costPlusZ - 2*cost + costMinZ ) / ( delta * delta );
  curvature.m_Curvature_dxdy = ( ( costPlusXPlusY - costPlusY ) - ( costPlusX - cost ) ) / ( delta * delta );
  curvature.m_Curvature_dxdz = ( ( costPlusXPlusZ - costPlusZ ) - ( costPlusX - cost ) ) / ( delta * delta );
  curvature.m_Curvature_dydz = ( ( costPlusYPlusZ - costPlusZ ) - ( costPlusY - cost ) ) / ( delta * delta );


  // Let's do a courtesy towards the user: if some components are immobile, the corresponding
  // rows and columns are now filled with zeros. Replace the elements on the diagonal with a
  // one, so that calculating the determinant yields the correct result
  if ( !m_MeshCollection->GetPointParameters()->ElementAt( pointId ).m_CanMoveX )
  {
    curvature.m_Curvature_dxdx = 1;
  }
  if ( !m_MeshCollection->GetPointParameters()->ElementAt( pointId ).m_CanMoveY )
  {
    curvature.m_Curvature_dydy = 1;
  }
  if ( !m_MeshCollection->GetPointParameters()->ElementAt( pointId ).m_CanMoveZ )
  {
    curvature.m_Curvature_dzdz = 1;
  }


  return curvature;

}


#endif


//
//
//
float
AtlasMeshVertexProcessor
::CalculatePriorCost( int meshNumber, AtlasMesh::PointIdentifier pointId, float x, float y, float z )
{
  // Make sure we precompute the vertex neighborhoods
  if ( m_VertexNeighborhoodsContainers.size() == 0 )
  {
    this->CalculateVertexNeighborhoods();
  }

  // Sanity checking
  if ( meshNumber >= m_VertexNeighborhoodsContainers.size() )
  {
    itkExceptionMacro( << "meshNumber out of range" );
  }
  if ( !( m_VertexNeighborhoodsContainers[ meshNumber ]->IndexExists( pointId ) ) )
  {
    itkExceptionMacro( << "pointId " << pointId << "does not exist" );
  }


  // Retrieve the neighborhood of this vertex
  const VertexNeighborhood&   neighborhood = m_VertexNeighborhoodsContainers[ meshNumber ]->ElementAt( pointId );

  // Loop over all tetrahedra, adding each tetrahedron's contribution to the cost
  float  cost = 0;
  for ( VertexNeighborhood::const_iterator tetrahedronInfoIt = neighborhood.begin();
        tetrahedronInfoIt != neighborhood.end();
        ++tetrahedronInfoIt )
  {
    // Get the information of the location of all four vertices, and fill them into
    // the matrix Y = [ p0 p1 p2 p3; 1 1 1 1 ];
    const float y11 = x;
    const float y21 = y;
    const float y31 = z;
    const float y12 = *( ( *tetrahedronInfoIt ).m_X1 );
    const float y22 = *( ( *tetrahedronInfoIt ).m_Y1 );
    const float y32 = *( ( *tetrahedronInfoIt ).m_Z1 );
    const float y13 = *( ( *tetrahedronInfoIt ).m_X2 );
    const float y23 = *( ( *tetrahedronInfoIt ).m_Y2 );
    const float y33 = *( ( *tetrahedronInfoIt ).m_Z2 );
    const float y14 = *( ( *tetrahedronInfoIt ).m_X3 );
    const float y24 = *( ( *tetrahedronInfoIt ).m_Y3 );
    const float y34 = *( ( *tetrahedronInfoIt ).m_Z3 );



    // Get the information we need about the reference position of this tetrahedron. It is
    // in the form of a precalculated Z = inv( Yref )
    const float  referenceVolumeTimesK = *( ( *tetrahedronInfoIt ).m_ReferenceVolumeTimesK );
    const float  z11 = *( ( *tetrahedronInfoIt ).m_Z11 );
    const float  z21 = *( ( *tetrahedronInfoIt ).m_Z21 );
    const float  z31 = *( ( *tetrahedronInfoIt ).m_Z31 );
    const float  z41 = *( ( *tetrahedronInfoIt ).m_Z41 );
    const float  z12 = *( ( *tetrahedronInfoIt ).m_Z12 );
    const float  z22 = *( ( *tetrahedronInfoIt ).m_Z22 );
    const float  z32 = *( ( *tetrahedronInfoIt ).m_Z32 );
    const float  z42 = *( ( *tetrahedronInfoIt ).m_Z42 );
    const float  z13 = *( ( *tetrahedronInfoIt ).m_Z13 );
    const float  z23 = *( ( *tetrahedronInfoIt ).m_Z23 );
    const float  z33 = *( ( *tetrahedronInfoIt ).m_Z33 );
    const float  z43 = *( ( *tetrahedronInfoIt ).m_Z43 );


    //
    // Now let's add Ashburner's prior cost for the tethrahedron deformation from its reference position
    //
    const float  m11 = z11*y11 + z21*y12 + z31*y13 + z41*y14;
    const float  m21 = z11*y21 + z21*y22 + z31*y23 + z41*y24;
    const float  m31 = z11*y31 + z21*y32 + z31*y33 + z41*y34;
    const float  m12 = z12*y11 + z22*y12 + z32*y13 + z42*y14;
    const float  m22 = z12*y21 + z22*y22 + z32*y23 + z42*y24;
    const float  m32 = z12*y31 + z22*y32 + z32*y33 + z42*y34;
    const float  m13 = z13*y11 + z23*y12 + z33*y13 + z43*y14;
    const float  m23 = z13*y21 + z23*y22 + z33*y23 + z43*y24;
    const float  m33 = z13*y31 + z23*y32 + z33*y33 + z43*y34;

    const float  detJ = m11 * ( m22*m33 - m32*m23 ) - m12 * ( m21*m33 - m31*m23 ) + m13 * ( m21*m32 - m31*m22 );
    if ( detJ <= 0 )
    {
      cost = itk::NumericTraits< float >::max();
      break;
    }


    // Let's define K as inv( J ) * det( J )
    const float  k11 = ( m22*m33 - m23*m32 );
    const float  k12 = -( m12*m33 - m32*m13 );
    const float  k13 = ( m12*m23 - m22*m13 );
    const float  k21 = -( m21*m33 - m31*m23 );
    const float  k22 = ( m11*m33 - m13*m31 );
    const float  k23 = -( m11*m23 - m21*m13 );
    const float  k31 = ( m21*m32 - m31*m22 );
    const float  k32 = -( m11*m32 - m31*m12 );
    const float  k33 = ( m11*m22 - m12*m21 );

    // Trace of J' * J is actually the sum of the squares of the singular values of J: s1^2 + s2^2 + s3^2
    const float  sumOfSquaresOfSingularValuesOfJ = m11*m11 + m12*m12 + m13*m13 + m21*m21 + m22*m22 + m23*m23 + m31*m31 + m32*m32 + m33*m33;

    // Trace of ( inv(J) )' * inv( J ) is actually the sum of the squares of the reciprocals of the singular values
    // of J: 1/s1^2 + 1/s2^2 + 1/s3^2
    const float  traceOfKTransposeTimesK = ( k11*k11 + k12*k12 + k13*k13 + k21*k21 + k22*k22 + k23*k23 + k31*k31 + k32*k32 + k33*k33 );
    const float  sumOfSquaresOfReciprocalsOfSingularValuesOfJ = traceOfKTransposeTimesK / ( detJ * detJ );

    const float  priorCost = referenceVolumeTimesK * ( 1 + detJ ) *
                             ( sumOfSquaresOfSingularValuesOfJ + sumOfSquaresOfReciprocalsOfSingularValuesOfJ - 6 );
    cost += priorCost;

  } // End loop over tetrahedra to calculate cost

  // Return result
  return cost;
}




//
//
//
AtlasPositionGradientType
AtlasMeshVertexProcessor
::CalculatePriorGradient( int meshNumber, AtlasMesh::PointIdentifier pointId, float x, float y, float z )
{
  // Make sure we precompute the vertex neighborhoods
  if ( m_VertexNeighborhoodsContainers.size() == 0 )
  {
    this->CalculateVertexNeighborhoods();
  }

  // Sanity checking
  if ( meshNumber >= m_VertexNeighborhoodsContainers.size() )
  {
    itkExceptionMacro( << "meshNumber out of range" );
  }
  if ( !( m_VertexNeighborhoodsContainers[ meshNumber ]->IndexExists( pointId ) ) )
  {
    itkExceptionMacro( << "pointId " << pointId << "does not exist" );
  }


  // Retrieve the neighborhood of this vertex
  const VertexNeighborhood&   neighborhood = m_VertexNeighborhoodsContainers[ meshNumber ]->ElementAt( pointId );


  // Calculate gradient by looping over all neighboring tetrahedra
  float dCostdx = 0;
  float dCostdy = 0;
  float dCostdz = 0;
  for ( VertexNeighborhood::const_iterator tetrahedronInfoIt = neighborhood.begin();
        tetrahedronInfoIt != neighborhood.end();
        ++tetrahedronInfoIt )
  {
    // Get the information of the location of all four vertices, and fill them into
    // the matrix Y = [ p0 p1 p2 p3; 1 1 1 1 ];
    const float y11 = x;
    const float y21 = y;
    const float y31 = z;
    const float y12 = *( ( *tetrahedronInfoIt ).m_X1 );
    const float y22 = *( ( *tetrahedronInfoIt ).m_Y1 );
    const float y32 = *( ( *tetrahedronInfoIt ).m_Z1 );
    const float y13 = *( ( *tetrahedronInfoIt ).m_X2 );
    const float y23 = *( ( *tetrahedronInfoIt ).m_Y2 );
    const float y33 = *( ( *tetrahedronInfoIt ).m_Z2 );
    const float y14 = *( ( *tetrahedronInfoIt ).m_X3 );
    const float y24 = *( ( *tetrahedronInfoIt ).m_Y3 );
    const float y34 = *( ( *tetrahedronInfoIt ).m_Z3 );


    // Get the information we need about the reference position of this tetrahedron. It is
    // in the form of a precalculated Z = inv( Yref )
    const float  referenceVolumeTimesK = *( ( *tetrahedronInfoIt ).m_ReferenceVolumeTimesK );
    const float  z11 = *( ( *tetrahedronInfoIt ).m_Z11 );
    const float  z21 = *( ( *tetrahedronInfoIt ).m_Z21 );
    const float  z31 = *( ( *tetrahedronInfoIt ).m_Z31 );
    const float  z41 = *( ( *tetrahedronInfoIt ).m_Z41 );
    const float  z12 = *( ( *tetrahedronInfoIt ).m_Z12 );
    const float  z22 = *( ( *tetrahedronInfoIt ).m_Z22 );
    const float  z32 = *( ( *tetrahedronInfoIt ).m_Z32 );
    const float  z42 = *( ( *tetrahedronInfoIt ).m_Z42 );
    const float  z13 = *( ( *tetrahedronInfoIt ).m_Z13 );
    const float  z23 = *( ( *tetrahedronInfoIt ).m_Z23 );
    const float  z33 = *( ( *tetrahedronInfoIt ).m_Z33 );
    const float  z43 = *( ( *tetrahedronInfoIt ).m_Z43 );


    //
    // Now let's add Ashburner's prior gradient for the tethrahedron deformation from its reference position
    //
    const float  m11 = z11*y11 + z21*y12 + z31*y13 + z41*y14;
    const float  m21 = z11*y21 + z21*y22 + z31*y23 + z41*y24;
    const float  m31 = z11*y31 + z21*y32 + z31*y33 + z41*y34;
    const float  m12 = z12*y11 + z22*y12 + z32*y13 + z42*y14;
    const float  m22 = z12*y21 + z22*y22 + z32*y23 + z42*y24;
    const float  m32 = z12*y31 + z22*y32 + z32*y33 + z42*y34;
    const float  m13 = z13*y11 + z23*y12 + z33*y13 + z43*y14;
    const float  m23 = z13*y21 + z23*y22 + z33*y23 + z43*y24;
    const float  m33 = z13*y31 + z23*y32 + z33*y33 + z43*y34;

    const float  detJ = m11 * ( m22*m33 - m32*m23 ) - m12 * ( m21*m33 - m31*m23 ) + m13 * ( m21*m32 - m31*m22 );

    // Let's define K as inv( J ) * det( J )
    const float  k11 = ( m22*m33 - m23*m32 );
    const float  k12 = -( m12*m33 - m32*m13 );
    const float  k13 = ( m12*m23 - m22*m13 );
    const float  k21 = -( m21*m33 - m31*m23 );
    const float  k22 = ( m11*m33 - m13*m31 );
    const float  k23 = -( m11*m23 - m21*m13 );
    const float  k31 = ( m21*m32 - m31*m22 );
    const float  k32 = -( m11*m32 - m31*m12 );
    const float  k33 = ( m11*m22 - m12*m21 );

    // Trace of J' * J is actually the sum of the squares of the singular values of J: s1^2 + s2^2 + s3^2
    const float  sumOfSquaresOfSingularValuesOfJ = m11*m11 + m12*m12 + m13*m13 + m21*m21 + m22*m22 + m23*m23 + m31*m31 + m32*m32 + m33*m33;

    // Trace of ( inv(J) )' * inv( J ) is actually the sum of the squares of the reciprocals of the singular values
    // of J: 1/s1^2 + 1/s2^2 + 1/s3^2
    const float  traceOfKTransposeTimesK = ( k11*k11 + k12*k12 + k13*k13 + k21*k21 + k22*k22 + k23*k23 + k31*k31 + k32*k32 + k33*k33 );
    const float  sumOfSquaresOfReciprocalsOfSingularValuesOfJ = traceOfKTransposeTimesK / ( detJ * detJ );


    //
    // OK, now add contribution to derivatives of Ashburner's prior cost in each of the tetrahedron's vertices
    //
    const float  ddetJdm11 = m22*m33 - m32*m23;
    const float  ddetJdm21 = m13*m32 - m12*m33;
    const float  ddetJdm31 = m12*m23 - m13*m22;
    const float  ddetJdm12 = m31*m23 - m21*m33;
    const float  ddetJdm22 = m11*m33 - m13*m31;
    const float  ddetJdm32 = m13*m21 - m11*m23;
    const float  ddetJdm13 = m21*m32 - m31*m22;
    const float  ddetJdm23 = m12*m31 - m11*m32;
    const float  ddetJdm33 = m11*m22 - m12*m21;


    const float  tmp1 = referenceVolumeTimesK * ( sumOfSquaresOfSingularValuesOfJ + sumOfSquaresOfReciprocalsOfSingularValuesOfJ - 6 );
    const float  tmp2 = 2 * referenceVolumeTimesK * ( 1 + detJ );
    const float  tmp3 = pow( detJ,  3 );


    const float  dcostdm11 = tmp1 * ddetJdm11 + tmp2 * ( m11 + ( ( k22*m33 - k23*m23 - k32*m32 + k33*m22 ) * detJ - traceOfKTransposeTimesK * ddetJdm11 ) / tmp3 );
    const float  dcostdm21 = tmp1 * ddetJdm21 + tmp2 * ( m21 + ( ( -k21*m33 + k23*m13 + k31*m32 - k33*m12 ) * detJ - traceOfKTransposeTimesK * ddetJdm21 ) / tmp3 );
    const float  dcostdm31 = tmp1 * ddetJdm31 + tmp2 * ( m31 + ( ( k21*m23 - k22*m13 - k31*m22 + k32*m12 ) * detJ - traceOfKTransposeTimesK * ddetJdm31 ) / tmp3 );

    const float  dcostdm12 = tmp1 * ddetJdm12 + tmp2 * ( m12 + ( ( -k12*m33 + k13*m23 + k32*m31 - k33*m21 ) * detJ - traceOfKTransposeTimesK * ddetJdm12 ) / tmp3 );
    const float  dcostdm22 = tmp1 * ddetJdm22 + tmp2 * ( m22 + ( ( k11*m33 - k13*m13 - k31*m31 + k33*m11 ) * detJ - traceOfKTransposeTimesK * ddetJdm22 ) / tmp3 );
    const float  dcostdm32 = tmp1 * ddetJdm32 + tmp2 * ( m32 + ( ( -k11*m23 + k12*m13 + k31*m21 - k32*m11 ) * detJ - traceOfKTransposeTimesK * ddetJdm32 ) / tmp3 );

    const float  dcostdm13 = tmp1 * ddetJdm13 + tmp2 * ( m13 + ( ( k12*m32 - k13*m22 - k22*m31 + k23*m21 ) * detJ - traceOfKTransposeTimesK * ddetJdm13 ) / tmp3 );
    const float  dcostdm23 = tmp1 * ddetJdm23 + tmp2 * ( m23 + ( ( -k11*m32 + k13*m12 + k21*m31 -k23*m11 ) * detJ - traceOfKTransposeTimesK * ddetJdm23 ) / tmp3 );
    const float  dcostdm33 = tmp1 * ddetJdm33 + tmp2 * ( m33 + ( ( k11*m22 - k12*m12 - k21*m21 + k22*m11 ) * detJ - traceOfKTransposeTimesK * ddetJdm33 ) / tmp3 );


    const float  dcostdy11 = dcostdm11*z11 + dcostdm12*z12 + dcostdm13*z13;
    const float  dcostdy21 = dcostdm21*z11 + dcostdm22*z12 + dcostdm23*z13;
    const float  dcostdy31 = dcostdm31*z11 + dcostdm32*z12 + dcostdm33*z13;

    // Add the stuff to the existing gradients in each of the tetrahedron's four vertices
    dCostdx += dcostdy11;
    dCostdy += dcostdy21;
    dCostdy += dcostdy31;

  } // End calculation of gradient by looping over all neighboring tetrahedra


  // Return result
  AtlasPositionGradientType  gradient;
  const bool  canMoveX = m_MeshCollection->GetPointParameters()->ElementAt( pointId ).m_CanMoveX;
  const bool  canMoveY = m_MeshCollection->GetPointParameters()->ElementAt( pointId ).m_CanMoveY;
  const bool  canMoveZ = m_MeshCollection->GetPointParameters()->ElementAt( pointId ).m_CanMoveZ;
  if ( !canMoveX )
  {
    gradient[ 0 ] = 0;
  }
  else
  {
    gradient[ 0 ] = dCostdx;
  }
  if ( !canMoveY )
  {
    gradient[ 1 ] = 0;
  }
  else
  {
    gradient[ 1 ] = dCostdy;
  }
  if ( !canMoveZ )
  {
    gradient[ 2 ] = 0;
  }
  else
  {
    gradient[ 2 ] = dCostdz;
  }

  return gradient;
}





#if 0

//
//
//
Curvature
AtlasMeshVertexProcessor
::CalculatePriorCurvature( int meshNumber, AtlasMesh::PointIdentifier pointId, float x, float y )
{
  // Make sure we precompute the vertex neighborhoods
  if ( m_VertexNeighborhoodsContainers.size() == 0 )
  {
    this->CalculateVertexNeighborhoods();
  }

  // Sanity checking
  if ( meshNumber >= m_VertexNeighborhoodsContainers.size() )
  {
    itkExceptionMacro( << "meshNumber out of range" );
  }
  if ( !( m_VertexNeighborhoodsContainers[ meshNumber ]->IndexExists( pointId ) ) )
  {
    itkExceptionMacro( << "pointId " << pointId << "does not exist" );
  }


  // Retrieve the neighborhood of this vertex
  const VertexNeighborhood&   neighborhood = m_VertexNeighborhoodsContainers[ meshNumber ]->ElementAt( pointId );

  // Calculate second derivates by looping over all neighboring triangles
  float d2Costdx2 = 0;
  float d2Costdxdy = 0;
  float d2Costdy2 = 0;
  for ( VertexNeighborhood::const_iterator triangleInfoIt = neighborhood.begin();
        triangleInfoIt != neighborhood.end();
        ++triangleInfoIt )
  {
    const float x1 = *( ( *triangleInfoIt ).m_X1 );
    const float y1 = *( ( *triangleInfoIt ).m_Y1 );
    const float x2 = *( ( *triangleInfoIt ).m_X2 );
    const float y2 = *( ( *triangleInfoIt ).m_Y2 );

    const float referenceAreaTimesK = *( ( *triangleInfoIt ).m_ReferenceAreaTimesK );
    const float z11 = *( ( *triangleInfoIt ).m_Z11 );
    const float z21 = *( ( *triangleInfoIt ).m_Z21 );
    const float z31 = *( ( *triangleInfoIt ).m_Z31 );
    const float z12 = *( ( *triangleInfoIt ).m_Z12 );
    const float z22 = *( ( *triangleInfoIt ).m_Z22 );
    const float z32 = *( ( *triangleInfoIt ).m_Z32 );

    // Calculate deformation field matrix M components
    const float  m11 = z11 * x + z21 * x1 + z31 * x2;
    const float  m12 = z12 * x + z22 * x1 + z32 * x2;
    const float  m21 = z11 * y + z21 * y1 + z31 * y2;
    const float  m22 = z12 * y + z22 * y1 + z32 * y2;

    // Precalculate some components that re-occur over and over again in the prior equations
    const float  detM = m11 * m22 - m12 * m21;
    const float  onePlusDetM = 1 + detM;
    const float  normM = m11*m11 + m12*m12 + m21*m21 + m22*m22;
    const float  oneOverDet3 = 1 / pow( detM, 3 );
    const float  oneOverDet4 = 1 / pow( detM, 4 );
    const float  tmp1OfM = 1 + 1 / pow( detM, 2 );
    const float  tmp2OfM = normM * tmp1OfM - 4;

    const float  dxFactor = 2 * ( m11 * z11 + m12*z12 ) * tmp1OfM -
                            2 * normM * ( z11 * m22 - z12 * m21 ) * oneOverDet3;
    const float  dyFactor = 2 * ( m21 * z11 + m22 * z12 ) * tmp1OfM -
                            2 * normM * ( m11 * z12 - m12 * z11 ) * oneOverDet3;

    d2Costdx2 += referenceAreaTimesK *
                 ( 2 * ( z11 * m22 - z12 * m21 ) * dxFactor +
                   onePlusDetM *
                   ( 2 * ( z11 * z11 + z12 * z12 ) * tmp1OfM -
                     8 * ( m11 * z11 + m12 * z12 ) * oneOverDet3 * ( z11 * m22 - z12 * m21 ) +
                     6 * normM * pow( z11 * m22 - z12 * m21, 2 ) * oneOverDet4 ) );

    d2Costdxdy += referenceAreaTimesK *
                  ( ( z11 * m22 - z12 * m21 ) * dyFactor +
                    ( z12 * m11 - z11 * m12 ) * dxFactor +
                    onePlusDetM *
                    ( (-4) * ( m11 * z11 + m12 * z12 ) * oneOverDet3 * ( m11 * z12 - m12 * z11 ) -
                      4  * ( m21 * z11 + m22 * z12 ) * oneOverDet3 * ( m22 * z11 - m21 * z12 ) +
                      6  * normM * ( m22 * z11 - m21 * z12 ) * ( m11 * z12 - m12 * z11 ) *
                      oneOverDet4 ) );

    d2Costdy2 += referenceAreaTimesK *
                 ( 2 * ( m11 * z12 - m12 * z11 ) * dyFactor +
                   onePlusDetM *
                   ( 2 * ( z11 * z11 + z12 * z12 ) * tmp1OfM -
                     8 * ( m21 * z11 + m22 * z12 ) * oneOverDet3 * ( m11 * z12 - m12 * z11 ) +
                     6 * normM * pow( m11 * z12 - m12 * z11, 2 ) * oneOverDet4 ) );

  } // End calculation of second derivatives by looping over all neighboring triangles


  // Return result
  Curvature  curvature;
  curvature.m_Curvature_dxdx = d2Costdx2;
  curvature.m_Curvature_dxdy = d2Costdxdy;
  curvature.m_Curvature_dydy = d2Costdy2;
  const bool  canMoveX = m_MeshCollection->GetPointParameters()->ElementAt( pointId ).m_CanMoveX;
  const bool  canMoveY = m_MeshCollection->GetPointParameters()->ElementAt( pointId ).m_CanMoveY;
  if ( !canMoveX )
  {
    curvature.m_Curvature_dxdx = itk::NumericTraits< float >::max();
    curvature.m_Curvature_dxdy = 0;
  }
  if ( !canMoveY )
  {
    curvature.m_Curvature_dydy = itk::NumericTraits< float >::max();
    curvature.m_Curvature_dxdy = 0;
  }

  return curvature;
}

#endif



//
//
//
bool
AtlasMeshVertexProcessor
::CalculateXstar( int meshNumber, AtlasMesh::PointIdentifier pointId,
                  float& xstar, float& ystar, float& zstar, bool verbose )
{

  // Make sure we precompute the vertex neighborhoods
  if ( m_VertexNeighborhoodsContainers.size() == 0 )
  {
    this->CalculateVertexNeighborhoods();
  }

  // Sanity checking
  if ( meshNumber >= m_VertexNeighborhoodsContainers.size() )
  {
    if ( verbose )
    {
      std::cout << "meshNumber out of range" << std::endl;
    }
    return false;
  }
  if ( !( m_VertexNeighborhoodsContainers[ meshNumber ]->IndexExists( pointId ) ) )
  {
    if ( verbose )
    {
      std::cout << "pointId " << pointId << "does not exist" << std::endl;
    }
    return false;
  }

  // Retrieve movability info of this point, and its current position
  const bool  canMoveX = m_MeshCollection->GetPointParameters()->ElementAt( pointId ).m_CanMoveX;
  const bool  canMoveY = m_MeshCollection->GetPointParameters()->ElementAt( pointId ).m_CanMoveY;
  const bool  canMoveZ = m_MeshCollection->GetPointParameters()->ElementAt( pointId ).m_CanMoveZ;
  const AtlasMesh::PointType&  position = m_MeshCollection->GetPositions()[ meshNumber ]->ElementAt( pointId );
  if ( !canMoveX && !canMoveY && !canMoveZ )
  {
    xstar = position[ 0 ];
    ystar = position[ 1 ];
    zstar = position[ 2 ];
    return true;
  }


  // Retrieve the neighborhood of this vertex
  const VertexNeighborhood&   neighborhood = m_VertexNeighborhoodsContainers[ meshNumber ]->ElementAt( pointId );


  //
  // Part 1. Calculate ( xstar, ystar, zstar )
  //

  // Initialize ( xstar, ystar, zstar ) with the original position in the mesh
  xstar = position[ 0 ];
  ystar = position[ 1 ];
  zstar = position[ 2 ];

  if ( verbose )
  {
    std::cout << "Initial position of (xstar, ystar, zstar ) for pointId " << pointId << " in meshNumber " << meshNumber
              << " : (" << xstar << ", " << ystar << ", " << zstar << ")" << std::endl;
  }



  // Estimate ( xstar, ystar, zstar ) using a Levenberg-Marquardt algorithm, i.e. look for the
  // vector step given by
  //
  //    ( J^t * J + lambda * diag( J^t * J ) ) * step =  J^t * ( - cost )
  //
  // where J = gradient^t
  //
  // For the practical implementation, we'll make use of the fact that the components of the
  // gradient that are immobile are automatically set to zero in this->CalculateGradient().
  // All we need to do is to replace, in the matrix "J^t * J + lambda * diag( J^t * J )", the
  // diagonal elements corresponding those immobile components by "1" (they are automatically
  // set to zero). If we do that, simply inverting the matrix will yield the correct results
  // for the other components, and simply using the usual equation for "step" (see above) will
  // then have a zero component for the immobile component because J^t will have a zero there.
  //
  // NOTE: I'm sure numerical people will turn around in their grave if they see I'm explicitly
  // inverting the matrices, but it should work all right for this low-dimensional problem.


  float  lambda = 1.0f;
  bool  converged = false;
  for ( int iterationNumber = 0; ( iterationNumber < 30 ) && ( !converged ); iterationNumber++ )
  {
    // Get cost, gradient, and curvature at the current position
    const float  currentCost = this->CalculateCost( meshNumber, pointId, xstar, ystar, zstar );
    const AtlasPositionGradientType  gradient = this->CalculateGradient( meshNumber, pointId, xstar, ystar, zstar );

#if 1
    if ( iterationNumber == 0 )
    {
      AtlasMeshPositionGradientCalculator::Pointer  gradientCalculator = AtlasMeshPositionGradientCalculator::New();
      gradientCalculator->SetLabelImage( m_LabelImages[ meshNumber ] );
      gradientCalculator->Rasterize( m_MeshCollection->GetMesh( meshNumber ) );
      const float  minLogLikelihoodTimesPrior = gradientCalculator->GetMinLogLikelihoodTimesPrior();
      std::cout << "                 currentCost: we SHOULD get " << minLogLikelihoodTimesPrior << std::endl;
      std::cout << "                 gradient: we SHOULD get" << gradientCalculator->GetPositionGradient()->ElementAt( pointId ) << std::endl;
    }
#endif




#if 1
    // Simple gradient descent
    const float  gradientMagnitude = gradient.GetNorm();
    const float  stepSize = 0.5;
    xstar -= gradient[ 0 ] / gradientMagnitude * stepSize;
    ystar -= gradient[ 1 ] / gradientMagnitude * stepSize;
    zstar -= gradient[ 2 ] / gradientMagnitude * stepSize;
    if ( verbose )
    {
      std::cout << "                 currentCost: " << currentCost << std::endl;
      std::cout << "                 gradient: " << gradient << std::endl;
      std::cout << "      Iteration " << iterationNumber << " -> (" << xstar
                << ", " << ystar << ", " << zstar<< ")" << std::endl;
    }
#else
    // Adjust lambda iteratively to find a balance between gradient descent and Newton's method that
    // brings in a better position
    while ( true )
    {
      // Get a tentative position

      // Construct the default LHS
      vnl_matrix_fixed< float, 3, 3 >  LHS;
      for ( int i = 0; i < 3; i++ )
      {
        for ( int j = 0; j < 3; j++ )
        {
          LHS( i, j ) = gradient[ i ] * gradient[ j ];
        }
        LHS( i, i ) *= ( 1 + lambda );
      }

      // Apply our trick
      for ( int i = 0; i < 3; i++ )
      {
        if ( gradient[ i ] == 0 )
        {
          LHS( i, i ) = 1;
        }
      }

      // Solve the equation for the step vector
      vnl_vector_fixed< float, 3 >  Jt;
      for ( int i = 0; i < 3; i++ )
      {
        Jt( i ) = gradient[ i ];
      }
      vnl_vector_fixed< float, 3 >  step = vnl_inverse( LHS ) * Jt;
      const float  xtentative = xstar - currentCost * step[ 0 ];
      const float  ytentative = ystar - currentCost * step[ 1 ];
      const float  ztentative = zstar - currentCost * step[ 2 ];


      // Get the cost at the tentative cost
      const float  tentativeCost = this->CalculateCost( meshNumber, pointId, xtentative, ytentative, ztentative );

      // Depending on the cost difference between the current position and the tentative position, decide what to do
      if ( tentativeCost <= currentCost )
      {
        xstar = xtentative;
        ystar = ytentative;
        zstar = ztentative;
        if ( verbose )
        {
          std::cout << "      Iteration " << iterationNumber << " -> (" << xstar
                    << ", " << ystar << ", " << zstar<< ")" << std::endl;
          std::cout << "               lambda: " << lambda << std::endl;
          std::cout << "                 cost: " << tentativeCost << std::endl;
        }

        // Check if converged
        if ( ( fabs( tentativeCost - currentCost ) / fabs( tentativeCost + currentCost ) * 2 ) <= 1E-10 )
        {
          converged = true;
          if ( verbose )
          {
            std::cout << "      Convergence detected (currentCost: " << currentCost
                      << ", tentativeCost: " << tentativeCost << ")" << std::endl;
          }
        }

        lambda /= 10;
        break;
      }
      else
      {
        //std::cout << "currentCost: " << currentCost << std::endl;
        //std::cout << "tentativeCost: " << tentativeCost << std::endl;
        //std::cout << "            dCostdx: " << dCostdx << std::endl;
        //std::cout << "            dCostdy: " << dCostdy << std::endl;
        //std::cout << "            dCostdz: " << dCostdz << std::endl;
        lambda *= 10;
        if ( lambda > 1E30 )
        {
          return false;
        }
        //std::cout << "==> Increasing lambda to " << lambda << std::endl;
      }


    } // End adjusting lambda until we move in the correct direction
#endif

  } // End Levenberg-Marquardt iterations


  return true;
}


} // end namespace kvl
