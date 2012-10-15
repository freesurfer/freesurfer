/**
 * @file  kvlAtlasMeshVertexProcessor2.cxx
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
#include "kvlAtlasMeshVertexProcessor2.h"

#include "kvlAtlasMeshRasterizor.h"
#include "kvlAtlasMeshPositionGradientCalculator.h"
#include "vnl/vnl_inverse.h"
#include "vnl/vnl_matrix_fixed.h"


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
float
AtlasMeshVertexProcessor
::CalculateCostAndGradient( int meshNumber, AtlasMesh::PointIdentifier pointId, float x, float y, float z,
                            AtlasPositionGradientType& gradient ) const
{

  // Sanity check on input
  if ( meshNumber >= static_cast< int >( m_LabelImages.size() ) )
  {
    itkExceptionMacro( << "Cannot determine data contribution to cost because the necessary label image"
                       " is not set" );
  }


  // Get the cell id of the vertex under investigation
  AtlasMesh::CellIdentifier  cellId = 0;
  const std::set< AtlasMesh::CellIdentifier >&  cellNeighbors = m_MeshCollection->GetCellLinks()->ElementAt( pointId );
  for ( std::set< AtlasMesh::CellIdentifier >::const_iterator  neighborIt = cellNeighbors.begin();
        neighborIt != cellNeighbors.end(); ++neighborIt )
  {
    const AtlasMesh::CellType*  cell = m_MeshCollection->GetCells()->ElementAt( *neighborIt );

    if ( cell->GetType() == AtlasMesh::CellType::VERTEX_CELL )
    {
      cellId = *neighborIt;
      break;
    }
  }
  //std::cout << "cellId: " << cellId << std::endl;


  // Get the mini mesh collection
  AtlasMeshCollection::Pointer  miniCollection = m_MeshCollection->GetRegionGrown( cellId, 1 );
  if ( !miniCollection )
  {
    itkExceptionMacro( << "unable to retrieve a mini mesh around point with id " << pointId );
  }


  // Get the mesh with the point under investigation set to the correct coordinates
  AtlasMesh::Pointer  mesh = const_cast< AtlasMesh* >( miniCollection->GetMesh( meshNumber ) );
  mesh->GetPoints()->ElementAt( pointId )[ 0 ] = x;
  mesh->GetPoints()->ElementAt( pointId )[ 1 ] = y;
  mesh->GetPoints()->ElementAt( pointId )[ 2 ] = z;


  // Set up a gradient calculator. Rather than simply calling
  //
  //   gradientCalculator->Rasterize( m_MeshCollection->GetMesh( meshNumber ) );
  //
  // on it, we're gonna abuse it by doing what happens internally in two separate steps.
  // This will allow us to calculate the prior and the non-normalized posterior parts
  // separately
  AtlasMeshPositionGradientCalculator::Pointer  gradientCalculator = AtlasMeshPositionGradientCalculator::New();
  gradientCalculator->SetLabelImage( m_LabelImages[ meshNumber ] );
  //gradientCalculator->Rasterize( m_MeshCollection->GetMesh( meshNumber ) );
  gradientCalculator->GetFragmentProcessor().SetMesh( mesh );


  // Loop over all tetrahedra, and pretend to start a new tetrahedron without ever rasterizing it.
  // This will allow to collect the prior part of the cost
  for ( AtlasMesh::CellsContainer::ConstIterator  cellIt = mesh->GetCells()->Begin();
        cellIt != mesh->GetCells()->End(); ++cellIt )
  {
    const AtlasMesh::CellType*  cell = cellIt.Value();

    if( cell->GetType() != AtlasMesh::CellType::TETRAHEDRON_CELL )
    {
      continue;
    }

    // OK, found tetrahedron. Warn the FragmentProcessor
    if ( !gradientCalculator->GetFragmentProcessor().StartNewTetrahedron( cellIt.Index() ) )
    {
      // Something wrong.
      break;
    }

  } // End loop over all tetrahedra

  const float  priorCost = gradientCalculator->GetMinLogLikelihoodTimesPrior();
  const AtlasPositionGradientType  priorGradient = gradientCalculator->GetPositionGradient()->ElementAt( pointId );
  //std::cout << "priorCost: " << priorCost << std::endl;
  //std::cout << "priorGradient: " << priorGradient << std::endl;

  // If beta is zero, that's all there is to it
  if ( ( m_Beta == 0 ) || ( priorCost == itk::NumericTraits< float >::max() ) )
  {
    gradient = priorGradient;
    if ( !m_MeshCollection->GetPointParameters()->ElementAt( pointId ).m_CanMoveX )
    {
      gradient[ 0 ] = 0;
    }
    if ( !m_MeshCollection->GetPointParameters()->ElementAt( pointId ).m_CanMoveY )
    {
      gradient[ 1 ] = 0;
    }
    if ( !m_MeshCollection->GetPointParameters()->ElementAt( pointId ).m_CanMoveZ )
    {
      gradient[ 2 ] = 0;
    }

    return priorCost;
  }


  // Beta is not zero. Rasterize the tetrahedra without collecting the prior part again
  gradientCalculator->GetFragmentProcessor().SetMesh( mesh ); // Will reset the cost to zero
  gradientCalculator->GetFragmentProcessor().SetIgnoreDeformationPrior( true );
  for ( AtlasMesh::CellsContainer::ConstIterator  cellIt = mesh->GetCells()->Begin();
        cellIt != mesh->GetCells()->End(); ++cellIt )
  {
    const AtlasMesh::CellType*  cell = cellIt.Value();

    if( cell->GetType() != AtlasMesh::CellType::TETRAHEDRON_CELL )
    {
      continue;
    }

    // OK, found tetrahedron. Warn the FragmentProcessor
    gradientCalculator->GetFragmentProcessor().StartNewTetrahedron( cellIt.Index() );

    // Retrieve position of vertices
    AtlasMesh::CellType::PointIdConstIterator  pit = cell->PointIdsBegin();
    AtlasMesh::PointType  p0;
    AtlasMesh::PointType  p1;
    AtlasMesh::PointType  p2;
    AtlasMesh::PointType  p3;
    mesh->GetPoint( *pit, &p0 );
    ++pit;
    mesh->GetPoint( *pit, &p1 );
    ++pit;
    mesh->GetPoint( *pit, &p2 );
    ++pit;
    mesh->GetPoint( *pit, &p3 );

    // Rasterize
    gradientCalculator->RasterizeTetrahedron( p0, p1, p2, p3 );

  } // End loop over all tetrahedra

  const float  dataCost = m_Beta * gradientCalculator->GetMinLogLikelihoodTimesPrior();
  const AtlasPositionGradientType  dataGradient = m_Beta * gradientCalculator->GetPositionGradient()->ElementAt( pointId );
  //std::cout << "dataCost: " << dataCost << std::endl;
  //std::cout << "dataGradient: " << dataGradient << std::endl;

  // Return sum of two contributions
  gradient = priorGradient + dataGradient;
  if ( !m_MeshCollection->GetPointParameters()->ElementAt( pointId ).m_CanMoveX )
  {
    gradient[ 0 ] = 0;
  }
  if ( !m_MeshCollection->GetPointParameters()->ElementAt( pointId ).m_CanMoveY )
  {
    gradient[ 1 ] = 0;
  }
  if ( !m_MeshCollection->GetPointParameters()->ElementAt( pointId ).m_CanMoveZ )
  {
    gradient[ 2 ] = 0;
  }

  return dataCost + priorCost;

}






//
//
//
float
AtlasMeshVertexProcessor
::CalculateCost( int meshNumber, AtlasMesh::PointIdentifier pointId, float x, float y, float z ) const
{
  AtlasPositionGradientType  dummyGradient;
  return this->CalculateCostAndGradient( meshNumber, pointId, x, y, z, dummyGradient );
}



//
//
//
AtlasPositionGradientType
AtlasMeshVertexProcessor
::CalculateGradient( int meshNumber, AtlasMesh::PointIdentifier pointId, float x, float y, float z ) const
{
  AtlasPositionGradientType  gradient;
  this->CalculateCostAndGradient( meshNumber, pointId, x, y, z, gradient );
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
::CalculateCurvature( int meshNumber, AtlasMesh::PointIdentifier pointId, float x, float y, float z, bool verbose ) const
{
  // Sanity check on input
  if ( meshNumber >= static_cast< int >( m_LabelImages.size() ) )
  {
    itkExceptionMacro( << "Cannot determine data contribution to cost because the necessary label image"
                       " is not set" );
  }


  // Quick and dirty: use finite differences
  const float  delta = 0.01 /* 1.0 */;
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
  if ( !( std::isnan( cost ) || std::isinf( cost ) || ( cost == itk::NumericTraits< float >::max() ) ||
          std::isnan( costPlusX ) || std::isinf( costPlusX ) || ( costPlusX == itk::NumericTraits< float >::max() ) ||
          std::isnan( costMinX ) || std::isinf( costMinX ) || ( costMinX == itk::NumericTraits< float >::max() ) ||
          std::isnan( costPlusY ) || std::isinf( costPlusY ) || ( costPlusY == itk::NumericTraits< float >::max() ) ||
          std::isnan( costMinY ) || std::isinf( costMinY ) || ( costMinY == itk::NumericTraits< float >::max() ) ||
          std::isnan( costPlusZ ) || std::isinf( costPlusZ ) || ( costPlusZ == itk::NumericTraits< float >::max() ) ||
          std::isnan( costMinZ ) || std::isinf( costMinZ ) || ( costMinZ == itk::NumericTraits< float >::max() ) ||
          std::isnan( costPlusXPlusY ) || std::isinf( costPlusXPlusY ) || ( costPlusXPlusY == itk::NumericTraits< float >::max() ) ||
          std::isnan( costPlusXPlusZ ) || std::isinf( costPlusXPlusZ ) || ( costPlusXPlusZ == itk::NumericTraits< float >::max() ) ||
          std::isnan( costPlusYPlusZ ) || std::isinf( costPlusYPlusZ ) || ( costPlusYPlusZ == itk::NumericTraits< float >::max() ) ) )
  {
    curvature.m_Curvature_dxdx = ( costPlusX - 2*cost + costMinX ) / ( delta * delta );
    curvature.m_Curvature_dydy = ( costPlusY - 2*cost + costMinY ) / ( delta * delta );
    curvature.m_Curvature_dzdz = ( costPlusZ - 2*cost + costMinZ ) / ( delta * delta );
    curvature.m_Curvature_dxdy = ( ( costPlusXPlusY - costPlusY ) - ( costPlusX - cost ) ) / ( delta * delta );
    curvature.m_Curvature_dxdz = ( ( costPlusXPlusZ - costPlusZ ) - ( costPlusX - cost ) ) / ( delta * delta );
    curvature.m_Curvature_dydz = ( ( costPlusYPlusZ - costPlusZ ) - ( costPlusY - cost ) ) / ( delta * delta );
  }


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
                  float& xstar, float& ystar, float& zstar, bool verbose ) const
{

  // Sanity check on input
  if ( ( meshNumber >= static_cast< int >( m_LabelImages.size() ) ) ||
       ( !m_MeshCollection->GetReferencePosition()->IndexExists( pointId ) ) )
  {
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
  AtlasPositionGradientType  currentGradient;
  float  currentCost = this->CalculateCostAndGradient( meshNumber, pointId, xstar, ystar, zstar,
                       currentGradient );
  for ( int iterationNumber = 0; ( iterationNumber < 30 ) && ( !converged ); iterationNumber++ )
  {
    //std::cout << "New iteration" << std::endl;

#if 0
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




#if 0
    // Simple gradient descent
    const float  gradientMagnitude = gradient.GetNorm() + 1e-16;
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
      // std::cout << "Going for round with lambda: " << lambda << std::endl;

      // Construct the default LHS
      vnl_matrix_fixed< float, 3, 3 >  LHS;
      for ( int i = 0; i < 3; i++ )
      {
        for ( int j = 0; j < 3; j++ )
        {
          LHS( i, j ) = currentGradient[ i ] * currentGradient[ j ];
        }
        LHS( i, i ) *= ( 1 + lambda );
      }

      // Apply our trick
      for ( int i = 0; i < 3; i++ )
      {
        if ( currentGradient[ i ] == 0 )
        {
          LHS( i, i ) = 1;
        }
      }

      // Solve the equation for the step vector
      vnl_vector_fixed< float, 3 >  Jt;
      for ( int i = 0; i < 3; i++ )
      {
        Jt( i ) = currentGradient[ i ];
      }
      vnl_vector_fixed< float, 3 >  step = vnl_inverse( LHS ) * Jt;
      const float  xtentative = xstar - currentCost * step[ 0 ];
      const float  ytentative = ystar - currentCost * step[ 1 ];
      const float  ztentative = zstar - currentCost * step[ 2 ];


      // Get the cost at the tentative cost
      // std::cout << "Getting tentative cost and gradient for position (" << xtentative << ", "
      //                                                                   << ytentative << ", "
      //                                                                   << ztentative << ")" << std::endl;
      // std::cout << "  currentCost: " << currentCost << std::endl;
      // std::cout << "  currentGradient: " << currentGradient << std::endl;
      // std::cout << "  LHS: " << LHS << std::endl;
      // std::cout << "  Jt: " << Jt << std::endl;
      AtlasPositionGradientType  tentativeGradient;
      float  tentativeCost = itk::NumericTraits< float >::max();
      if ( !( std::isnan( xtentative ) || std::isnan( ytentative ) || std::isnan( ztentative ) ||
              std::isinf( xtentative ) || std::isinf( ytentative ) || std::isinf( ztentative ) ) )
      {
        tentativeCost = this->CalculateCostAndGradient( meshNumber, pointId, xtentative, ytentative, ztentative,
                        tentativeGradient );
      }
      // std::cout << "Still here" << std::endl;

      // Depending on the cost difference between the current position and the tentative position, decide what to do
      if ( tentativeCost <= currentCost )
      {
        // Tentative step was good. Update the current position to it.
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

        currentCost = tentativeCost;
        currentGradient = tentativeGradient;

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
        //std::cout << "==> Increasing lambda to " << lambda << std::endl;
        lambda *= 10;
        if ( lambda > 1E30 )
        {
          return false;
        }
      }


    } // End adjusting lambda until we move in the correct direction
#endif

  } // End Levenberg-Marquardt iterations


  return true;
}


} // end namespace kvl
