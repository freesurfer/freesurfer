/**
 * @file  kvlAtlasMeshCollectionReferencePositionCost.cxx
 * @brief REPLACE_WITH_ONE_LINE_SHORT_DESCRIPTION
 *
 * REPLACE_WITH_LONG_DESCRIPTION_OR_REFERENCE
 */
/*
 * Original Author: Koen Van Leemput
 * CVS Revision Info:
 *    $Author: nicks $
 *    $Date: 2012/10/15 21:17:38 $
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
#include "kvlAtlasMeshCollectionReferencePositionCost.h"



namespace kvl
{


//
//
//
AtlasMeshCollectionReferencePositionCost
::AtlasMeshCollectionReferencePositionCost()
{

  m_Estimator = AtlasParameterEstimator::New();
  m_DataAndAlphasCostCalculator = AtlasMeshCollectionModelLikelihoodCalculator::New();
  m_PositionCostCalculator = AtlasMeshCollectionPositionCostCalculator::New();

  m_PointId = 0;
  m_InitialPointParameters = 0;
  m_Cells = 0;
  m_ReferencePosition = 0;
  m_K = 1;

  m_CostCalculationMeshCollection = 0;
}


//
//
//
AtlasMeshCollectionReferencePositionCost
::~AtlasMeshCollectionReferencePositionCost()
{
}



//
//
//
void
AtlasMeshCollectionReferencePositionCost
::SetInitialMeshCollections( AtlasMeshCollection* estimationMeshCollection,
                             AtlasMeshCollection* costCalculationMeshCollection )
{

  //
  m_InitialPositions = estimationMeshCollection->GetPositions();
  m_InitialPointParameters = estimationMeshCollection->GetPointParameters();
  m_Cells = estimationMeshCollection->GetCells();
  m_ReferencePosition = estimationMeshCollection->GetReferencePosition();
  m_K = estimationMeshCollection->GetK();

  // Any existing vertex neighborhood is no longer valid
  m_VertexNeighborhood.clear();

  m_CostCalculationMeshCollection = costCalculationMeshCollection;
}



//
//
//
void
AtlasMeshCollectionReferencePositionCost
::SetLabelImages( const std::vector< LabelImageType::ConstPointer >& labelImages )
{
  m_Estimator->SetLabelImages( labelImages );
  m_DataAndAlphasCostCalculator->SetLabelImages( labelImages );
  m_PositionCostCalculator->SetLabelImages( labelImages );
}





//
//
//
void
AtlasMeshCollectionReferencePositionCost
::
CalculateVertexNeighborhood() const
{
  // Sanity check
  if ( !m_ReferencePosition || !m_Cells )
  {
    itkExceptionMacro( "No reference position!" );
  }

  // Clean up
  m_VertexNeighborhood.clear();

  // Loop over all tetrahedra
  AtlasMesh::CellsContainer::ConstIterator  cellIt = m_Cells->Begin();
  while ( cellIt != m_Cells->End() )
  {
    const AtlasMesh::CellType*  cell = cellIt.Value();

    if( cell->GetType() != AtlasMesh::CellType::TETRAHEDRON_CELL )
    {
      ++cellIt;
      continue;
    }

    // Retrieve point id's of the four vertices
    AtlasMesh::CellType::PointIdConstIterator pointIt = cell->PointIdsBegin();
    AtlasMesh::PointIdentifier point0Id = *pointIt;
    ++pointIt;
    AtlasMesh::PointIdentifier point1Id = *pointIt;
    ++pointIt;
    AtlasMesh::PointIdentifier point2Id = *pointIt;
    ++pointIt;
    AtlasMesh::PointIdentifier point3Id = *pointIt;

    // We're only interested in those tetrahedra that have m_PointId as one of its vertices
    if ( ( point0Id != m_PointId ) && ( point1Id != m_PointId ) && ( point2Id != m_PointId ) && ( point3Id != m_PointId ) )
    {
      ++cellIt;
      continue;
    }


    // Retrieve address of (x, y, z) values in each of the vertices
    const AtlasMesh::PointType::ValueType  x0 = m_ReferencePosition->ElementAt( point0Id ).GetElement( 0 );
    const AtlasMesh::PointType::ValueType  y0 = m_ReferencePosition->ElementAt( point0Id ).GetElement( 1 );
    const AtlasMesh::PointType::ValueType  z0 = m_ReferencePosition->ElementAt( point0Id ).GetElement( 2 );
    const AtlasMesh::PointType::ValueType  x1 = m_ReferencePosition->ElementAt( point1Id ).GetElement( 0 );
    const AtlasMesh::PointType::ValueType  y1 = m_ReferencePosition->ElementAt( point1Id ).GetElement( 1 );
    const AtlasMesh::PointType::ValueType  z1 = m_ReferencePosition->ElementAt( point1Id ).GetElement( 2 );
    const AtlasMesh::PointType::ValueType  x2 = m_ReferencePosition->ElementAt( point2Id ).GetElement( 0 );
    const AtlasMesh::PointType::ValueType  y2 = m_ReferencePosition->ElementAt( point2Id ).GetElement( 1 );
    const AtlasMesh::PointType::ValueType  z2 = m_ReferencePosition->ElementAt( point2Id ).GetElement( 2 );
    const AtlasMesh::PointType::ValueType  x3 = m_ReferencePosition->ElementAt( point3Id ).GetElement( 0 );
    const AtlasMesh::PointType::ValueType  y3 = m_ReferencePosition->ElementAt( point3Id ).GetElement( 1 );
    const AtlasMesh::PointType::ValueType  z3 = m_ReferencePosition->ElementAt( point3Id ).GetElement( 2 );

    if ( point0Id == m_PointId )
    {
      // Add this tetrahedron to the neighborhood of vertex 0
      VertexNeighboringTetrahedronInfo  vertex0Info;
      vertex0Info.m_TetrahedronId = cellIt.Index();
      vertex0Info.m_X1 = x1;
      vertex0Info.m_Y1 = y1;
      vertex0Info.m_Z1 = z1;
      vertex0Info.m_X2 = x2;
      vertex0Info.m_Y2 = y2;
      vertex0Info.m_Z2 = z2;
      vertex0Info.m_X3 = x3;
      vertex0Info.m_Y3 = y3;
      vertex0Info.m_Z3 = z3;
      m_VertexNeighborhood.push_back( vertex0Info );
    }
    else if ( point1Id == m_PointId )
    {
      // Add this tetrahedron to the neighborhood of vertex 1
      VertexNeighboringTetrahedronInfo  vertex1Info;
      vertex1Info.m_TetrahedronId = cellIt.Index();
      vertex1Info.m_X1 = x2;
      vertex1Info.m_Y1 = y2;
      vertex1Info.m_Z1 = z2;
      vertex1Info.m_X2 = x0;
      vertex1Info.m_Y2 = y0;
      vertex1Info.m_Z2 = z0;
      vertex1Info.m_X3 = x3;
      vertex1Info.m_Y3 = y3;
      vertex1Info.m_Z3 = z3;
      m_VertexNeighborhood.push_back( vertex1Info );
    }
    else if ( point2Id == m_PointId )
    {
      // Add this tetrahedron to the neighborhood of vertex 2
      VertexNeighboringTetrahedronInfo  vertex2Info;
      vertex2Info.m_TetrahedronId = cellIt.Index();
      vertex2Info.m_X1 = x0;
      vertex2Info.m_Y1 = y0;
      vertex2Info.m_Z1 = z0;
      vertex2Info.m_X2 = x1;
      vertex2Info.m_Y2 = y1;
      vertex2Info.m_Z2 = z1;
      vertex2Info.m_X3 = x3;
      vertex2Info.m_Y3 = y3;
      vertex2Info.m_Z3 = z3;
      m_VertexNeighborhood.push_back( vertex2Info );
    }
    else
    {
      // Add this tetrahedron to the neighborhood of vertex 3
      VertexNeighboringTetrahedronInfo  vertex3Info;
      vertex3Info.m_TetrahedronId = cellIt.Index();
      vertex3Info.m_X1 = x2;
      vertex3Info.m_Y1 = y2;
      vertex3Info.m_Z1 = z2;
      vertex3Info.m_X2 = x1;
      vertex3Info.m_Y2 = y1;
      vertex3Info.m_Z2 = z1;
      vertex3Info.m_X3 = x0;
      vertex3Info.m_Y3 = y0;
      vertex3Info.m_Z3 = z0;
      m_VertexNeighborhood.push_back( vertex3Info );
    }


    ++cellIt;
  } // End loop over all tetrahedra




}


//
//
//
AtlasMeshCollectionReferencePositionCost::MeasureType
AtlasMeshCollectionReferencePositionCost
::GetValue( const ParametersType & parameters ) const
{
  // Get the individual cost components
  float  dataCost;
  float  alphasCost;
  float  positionCost;
  if ( !this->GetValue( parameters, dataCost, alphasCost, positionCost ) )
  {
    return itk::NumericTraits< float >::max();
  }

  // Add the components
  float  totalCost = dataCost + alphasCost + positionCost;
  if ( std::isnan( totalCost ) )
  {
    totalCost = itk::NumericTraits< float >::max();
  }

  // Return the result
  return totalCost;

}



//
//
//
bool
AtlasMeshCollectionReferencePositionCost
::GetValue( const ParametersType& parameters, float& dataCost, float& alphasCost, float& positionCost ) const
{

  // Make sure the vertex neighborhood has been calculated
  if ( m_VertexNeighborhood.size() == 0 )
  {
    this->CalculateVertexNeighborhood();
  }

  // Set the reference position in the vertex under study to the provided parameters
  m_ReferencePosition->ElementAt( m_PointId )[ 0 ] = parameters[ 0 ];
  m_ReferencePosition->ElementAt( m_PointId )[ 1 ] = parameters[ 1 ];
  m_ReferencePosition->ElementAt( m_PointId )[ 2 ] = parameters[ 2 ];


  //
  // Part I: analyze reference position w.r.t. rest of reference position
  //
  const float x0 = parameters[ 0 ];
  const float y0 = parameters[ 1 ];
  const float z0 = parameters[ 2 ];
  for ( std::vector< VertexNeighboringTetrahedronInfo >::const_iterator tetrahedronInfoIt = m_VertexNeighborhood.begin();
        tetrahedronInfoIt != m_VertexNeighborhood.end();
        ++tetrahedronInfoIt )
  {
    const float x1 = ( *tetrahedronInfoIt ).m_X1;
    const float y1 = ( *tetrahedronInfoIt ).m_Y1;
    const float z1 = ( *tetrahedronInfoIt ).m_Z1;
    const float x2 = ( *tetrahedronInfoIt ).m_X2;
    const float y2 = ( *tetrahedronInfoIt ).m_Y2;
    const float z2 = ( *tetrahedronInfoIt ).m_Z2;
    const float x3 = ( *tetrahedronInfoIt ).m_X3;
    const float y3 = ( *tetrahedronInfoIt ).m_Y3;
    const float z3 = ( *tetrahedronInfoIt ).m_Z3;



    // Test if this tethrahedron is left-rotating. If it's not, stop looping over all affected tetrahedra.
    // Do this by calculating the volume of the tetrahedron; this should be positive.
    // In what follows, the matrix Lambda is the Jacobian of the transform from a standarized tetrahedron
    // ( ( 0 0 0 )^T, ( 1 0 0 )^T, ( 0 1 0 )^T, ( 0 0 1 )^T ), which has volume 1/6, to the actual tetrahedron
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
      std::cout << "referencePosition would cause a tetrahedron flip; rejecting this possibility" << std::endl;
      return false;
    }
#if 1
    {
      // Do an Ashuburner-like cost for the reference position from the standard tetrahedron. We just
      // leave out the absolute volume component; only deformation beyond scaling is checked

      // Let's define K as inv( Lambda ) * det( Lambda )
      const float  k11 = ( lambda22*lambda33 - lambda23*lambda32 );
      const float  k12 = -( lambda12*lambda33 - lambda32*lambda13 );
      const float  k13 = ( lambda12*lambda23 - lambda22*lambda13 );
      const float  k21 = -( lambda21*lambda33 - lambda31*lambda23 );
      const float  k22 = ( lambda11*lambda33 - lambda13*lambda31 );
      const float  k23 = -( lambda11*lambda23 - lambda21*lambda13 );
      const float  k31 = ( lambda21*lambda32 - lambda31*lambda22 );
      const float  k32 = -( lambda11*lambda32 - lambda31*lambda12 );
      const float  k33 = ( lambda11*lambda22 - lambda12*lambda21 );

      // Trace of Lambda' * Lambda is actually the sum of the squares of the singular values of Lambda: s1^2 + s2^2 + s3^2
      //const float  sumOfSquaresOfSingularValuesOfLambda =
      //     lambda11^2 + lambda12^2 + lambda13^2 + lambda21^2 + lambda22^2 + lambda23^2 + lambda31^2 + lambda32^2 + lambda33^2;

      // Trace of ( inv(Lambda) )' * inv( Lambda ) is actually the sum of the squares of the reciprocals of the singular
      // values of Lambda: 1/s1^2 + 1/s2^2 + 1/s3^2
      const float  traceOfKTransposeTimesK = ( k11*k11 + k12*k12 + k13*k13 +
                                             k21*k21 + k22*k22 + k23*k23 +
                                             k31*k31 + k32*k32 + k33*k33 );
      const float  sumOfSquaresOfReciprocalsOfSingularValuesOfLambda = traceOfKTransposeTimesK / pow( 6 * volume, 2 );

      // badness is ( (s1*s2*s3)^(1/3) / s1 )^2 + ( (s1*s2*s3)^(1/3) / s2 )^2 + ( (s1*s2*s3)^(1/3) / s3 )^2 - 3
      // This is insensitive to scaling
      const float  badness = pow( 6*volume, 2.0f / 3.0f ) * sumOfSquaresOfReciprocalsOfSingularValuesOfLambda - 3;
      //std::cout << "badness is " << badness << std::endl;
      if ( badness > /*20 */ 100000000 )
      {
        std::cout << "referencePosition is too bad; rejecting this possibility" << std::endl;
        return false;
      }

    }
#endif



  }


  //
  // Part II: Construct a mesh collection, and estimate it.
  //

  // Initialize the positions to the initial positions
  std::vector< AtlasMesh::PointsContainer::Pointer >   positions;
  for ( unsigned int i = 0; i < m_InitialPositions.size(); i++ )
  {
    AtlasMesh::PointsContainer::Pointer  position = AtlasMesh::PointsContainer::New();
    AtlasMesh::PointsContainer::ConstPointer  initialPosition = m_InitialPositions[ i ].GetPointer();
    for ( AtlasMesh::PointsContainer::ConstIterator  initialIt = initialPosition->Begin();
          initialIt != initialPosition->End(); ++initialIt )
    {
      position->InsertElement( initialIt.Index(), initialIt.Value() );
    }

    positions.push_back( position );
  }

  // Initialize the alphas to the initial alphas
  AtlasMesh::PointDataContainer::Pointer  pointParameters = AtlasMesh::PointDataContainer::New();
  for ( AtlasMesh::PointDataContainer::ConstIterator  initialIt = m_InitialPointParameters->Begin();
        initialIt != m_InitialPointParameters->End(); ++initialIt )
  {
    pointParameters->InsertElement( initialIt.Index(), initialIt.Value() );
  }

  // Create a mesh collection
  AtlasMeshCollection::Pointer  meshCollection = AtlasMeshCollection::New();
  meshCollection->SetPositions( positions );
  meshCollection->SetPointParameters( pointParameters );
  meshCollection->SetCells( m_Cells );
  meshCollection->SetReferencePosition( m_ReferencePosition );
  meshCollection->SetK( m_K );

  // Estimate the mesh collection alphas and positions
  m_Estimator->SetInitialMeshCollection( meshCollection );
  m_Estimator->Estimate();
  //meshCollection->Write( "tmp.txt" );


  //
  // Part III. Copy the changes into the costCalculationMeshCollection, and calculate the cost
  //

  // Copy the reference position
  m_CostCalculationMeshCollection->GetReferencePosition()->ElementAt( m_PointId )[ 0 ] = parameters[ 0 ];
  m_CostCalculationMeshCollection->GetReferencePosition()->ElementAt( m_PointId )[ 1 ] = parameters[ 1 ];
  m_CostCalculationMeshCollection->GetReferencePosition()->ElementAt( m_PointId )[ 2 ] = parameters[ 2 ];

  // Make sure the triangle areas are updated!
  m_CostCalculationMeshCollection->SetK( m_CostCalculationMeshCollection->GetK() );

  // Copy the positions
  for ( unsigned int i = 0; i < positions.size(); i++ )
  {
    AtlasMesh::PointsContainer::ConstPointer  source = positions[ i ].GetPointer();
    AtlasMesh::PointsContainer::Pointer  target = m_CostCalculationMeshCollection->GetPositions()[ i ];
    for ( AtlasMesh::PointsContainer::ConstIterator  sourceIt = source->Begin();
          sourceIt != source->End(); ++sourceIt )
    {
      target->ElementAt( sourceIt.Index() ) = sourceIt.Value();
    }

  } // End loop over all positions

  // Copy the alphas
  for ( AtlasMesh::PointDataContainer::ConstIterator  sourceIt = pointParameters->Begin();
        sourceIt != pointParameters->End(); ++sourceIt )
  {
    m_CostCalculationMeshCollection->GetPointParameters()->ElementAt( sourceIt.Index() ).m_Alphas =
      sourceIt.Value().m_Alphas;
  }

  // Compute the total model cost
  m_DataAndAlphasCostCalculator->SetMeshCollection( m_CostCalculationMeshCollection );
  m_DataAndAlphasCostCalculator->GetDataCostAndAlphasCost( dataCost, alphasCost );
  m_PositionCostCalculator->SetMeshCollection( m_CostCalculationMeshCollection );
  positionCost = m_PositionCostCalculator->GetPositionCost();

#if 0
  if ( ( ( parameters[ 0 ] == 3 ) && ( parameters[ 1 ] == 9 ) ) ||
       ( ( parameters[ 0 ] == 4 ) && ( parameters[ 1 ] == 10 ) ) )
  {
    std::cout << "m_CostCalculationMeshCollection is now for " << parameters[0] << " " << parameters[1] << std::endl;
    std::ostringstream  fileNameStream;
    fileNameStream << "m_CostCalculationMeshCollection_" << parameters[0] << "_" << parameters[1] << ".txt";
    m_CostCalculationMeshCollection->Write( fileNameStream.str().c_str() );
  }
#endif

  return true;
}




} // End namespace kvl
