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
  m_mapCompToComp = 0;
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
  //// std::cout << "cellId: " << cellId << std::endl;


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
  gradientCalculator->SetMapCompToComp( m_mapCompToComp );
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
  //// std::cout << "priorCost: " << priorCost << std::endl;
  //// std::cout << "priorGradient: " << priorGradient << std::endl;

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
  gradientCalculator->SetMapCompToComp( m_mapCompToComp ); // I don't think I need this, but in case...
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
  const AtlasPositionGradientType  dataGradient = static_cast< double >( m_Beta ) * gradientCalculator->GetPositionGradient()->ElementAt( pointId );
  //// std::cout << "dataCost: " << dataCost << std::endl;
  //// std::cout << "dataGradient: " << dataGradient << std::endl;

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
    // std::cout << "cost: " << cost << std::endl;
    // std::cout << "costPlusX - cost: " << costPlusX - cost << std::endl;
    // std::cout << "costMinX - cost: " << costMinX - cost << std::endl;
    // std::cout << "costPlusY - cost: " << costPlusY - cost << std::endl;
    // std::cout << "costMinY - cost: " << costMinY - cost << std::endl;
    // std::cout << "costPlusZ - cost: " << costPlusZ - cost << std::endl;
    // std::cout << "costMinZ - cost: " << costMinZ - cost << std::endl;
    // std::cout << "costPlusXPlusY - cost: " << costPlusXPlusY - cost << std::endl;
    // std::cout << "costPlusXPlusZ - cost: " << costPlusXPlusZ - cost << std::endl;
    // std::cout << "costPlusYPlusZ - cost: " << costPlusYPlusZ - cost << std::endl;
    }


  Curvature  curvature;
  if ( !( isnan( cost ) || isinf( cost ) || ( cost == itk::NumericTraits< float >::max() ) ||
          isnan( costPlusX ) || isinf( costPlusX ) || ( costPlusX == itk::NumericTraits< float >::max() ) ||
          isnan( costMinX ) || isinf( costMinX ) || ( costMinX == itk::NumericTraits< float >::max() ) ||
          isnan( costPlusY ) || isinf( costPlusY ) || ( costPlusY == itk::NumericTraits< float >::max() ) ||
          isnan( costMinY ) || isinf( costMinY ) || ( costMinY == itk::NumericTraits< float >::max() ) ||
          isnan( costPlusZ ) || isinf( costPlusZ ) || ( costPlusZ == itk::NumericTraits< float >::max() ) ||
          isnan( costMinZ ) || isinf( costMinZ ) || ( costMinZ == itk::NumericTraits< float >::max() ) ||
          isnan( costPlusXPlusY ) || isinf( costPlusXPlusY ) || ( costPlusXPlusY == itk::NumericTraits< float >::max() ) ||
          isnan( costPlusXPlusZ ) || isinf( costPlusXPlusZ ) || ( costPlusXPlusZ == itk::NumericTraits< float >::max() ) ||
          isnan( costPlusYPlusZ ) || isinf( costPlusYPlusZ ) || ( costPlusYPlusZ == itk::NumericTraits< float >::max() ) ) )
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
    // std::cout << "Initial position of (xstar, ystar, zstar ) for pointId " << pointId << " in meshNumber " << meshNumber
    //          << " : (" << xstar << ", " << ystar << ", " << zstar << ")" << std::endl;
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
    //// std::cout << "New iteration" << std::endl;


#if 0
    // Simple gradient descent
    const float  gradientMagnitude = gradient.GetNorm() + 1e-16;
    const float  stepSize = 0.5;
    xstar -= gradient[ 0 ] / gradientMagnitude * stepSize;
    ystar -= gradient[ 1 ] / gradientMagnitude * stepSize;
    zstar -= gradient[ 2 ] / gradientMagnitude * stepSize;
    if ( verbose )
      {
      // std::cout << "                 currentCost: " << currentCost << std::endl;
      // std::cout << "                 gradient: " << gradient << std::endl;
      // std::cout << "      Iteration " << iterationNumber << " -> (" << xstar
                << ", " << ystar << ", " << zstar<< ")" << std::endl;
      }
#else
    // Adjust lambda iteratively to find a balance between gradient descent and Newton's method that
    // brings in a better position
    while ( true )
      {
      // Get a tentative position
      // // std::cout << "Going for round with lambda: " << lambda << std::endl;

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
      // // std::cout << "Getting tentative cost and gradient for position (" << xtentative << ", "
      //                                                                   << ytentative << ", "
      //                                                                   << ztentative << ")" << std::endl;
      // // std::cout << "  currentCost: " << currentCost << std::endl;
      // // std::cout << "  currentGradient: " << currentGradient << std::endl;
      // // std::cout << "  LHS: " << LHS << std::endl;
      // // std::cout << "  Jt: " << Jt << std::endl;
      AtlasPositionGradientType  tentativeGradient;
      float  tentativeCost = itk::NumericTraits< float >::max();
      if ( !( isnan( xtentative ) || isnan( ytentative ) || isnan( ztentative ) ||
              isinf( xtentative ) || isinf( ytentative ) || isinf( ztentative ) ) )
        {
        tentativeCost = this->CalculateCostAndGradient( meshNumber, pointId, xtentative, ytentative, ztentative,
                                                                   tentativeGradient );
        }
      // // std::cout << "Still here" << std::endl;

      // Depending on the cost difference between the current position and the tentative position, decide what to do
      if ( tentativeCost <= currentCost )
        {
        // Tentative step was good. Update the current position to it.
        xstar = xtentative;
        ystar = ytentative;
        zstar = ztentative;

        if ( verbose )
          {
          // std::cout << "      Iteration " << iterationNumber << " -> (" << xstar
          //                                                   << ", " << ystar << ", " << zstar<< ")" << std::endl;
          // std::cout << "               lambda: " << lambda << std::endl;
          // std::cout << "                 cost: " << tentativeCost << std::endl;
          }

        // Check if converged 
        if ( ( fabs( tentativeCost - currentCost ) / fabs( tentativeCost + currentCost ) * 2 ) <= 1E-10 )
          {
          converged = true;
          if ( verbose )
            {
            // std::cout << "      Convergence detected (currentCost: " << currentCost 
             //         << ", tentativeCost: " << tentativeCost << ")" << std::endl;
            }
          }

        currentCost = tentativeCost;
        currentGradient = tentativeGradient;

        lambda /= 10;
        break;
        }
      else
        {
        //// std::cout << "currentCost: " << currentCost << std::endl;
        //// std::cout << "tentativeCost: " << tentativeCost << std::endl;
        //// std::cout << "            dCostdx: " << dCostdx << std::endl;
        //// std::cout << "            dCostdy: " << dCostdy << std::endl;
        //// std::cout << "            dCostdz: " << dCostdz << std::endl;
        //// std::cout << "==> Increasing lambda to " << lambda << std::endl;
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
