#include "kvlAtlasMeshCollectionFastReferencePositionCost.h"

#include "kvlTetrahedronAspectRatio.h"


namespace kvl
{


//
//
//
AtlasMeshCollectionFastReferencePositionCost
::AtlasMeshCollectionFastReferencePositionCost()
{

  m_Estimator = AtlasParameterEstimator::New();
  m_Estimator->SetPositionOptimizer( AtlasParameterEstimator::FIXED_STEP_GRADIENT_DESCENT );
  m_Estimator->SetMaximumNumberOfIterations( 3 );
  m_Estimator->SetNumberOfThreads( 1 );
  //  m_Estimator->GetDeformationOptimizer()->SetVerbose( false );
  
  m_DataAndAlphasCostCalculator = AtlasMeshCollectionModelLikelihoodCalculator::New();

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
AtlasMeshCollectionFastReferencePositionCost
::~AtlasMeshCollectionFastReferencePositionCost()
{
}



//
//
//
void
AtlasMeshCollectionFastReferencePositionCost
::SetInitialMeshCollection( AtlasMeshCollection* meshCollection )
{

  // 
  m_InitialPositions = meshCollection->GetPositions();
  m_InitialPointParameters = meshCollection->GetPointParameters();
  m_Cells = meshCollection->GetCells();
  m_ReferencePosition = meshCollection->GetReferencePosition();
  m_K = meshCollection->GetK();
  
  // Any existing vertex neighborhood is no longer valid
  m_VertexNeighborhood.clear();

  m_CostCalculationMeshCollection = 0;
}



//
//
//
void
AtlasMeshCollectionFastReferencePositionCost
::SetLabelImages( const std::vector< LabelImageType::ConstPointer >& labelImages,
                  const CompressionLookupTable*  compressionLookupTable )
{
  m_Estimator->SetLabelImages( labelImages, compressionLookupTable );
  m_DataAndAlphasCostCalculator->SetLabelImages( labelImages, compressionLookupTable );

}
 




//
//
//
void
AtlasMeshCollectionFastReferencePositionCost
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
AtlasMeshCollectionFastReferencePositionCost::MeasureType
AtlasMeshCollectionFastReferencePositionCost
::GetValue( const ParametersType & parameters ) const
{
  // Get the individual cost components
  double  dataCost;
  double  alphasCost;
  double  positionCost;
  if ( !this->GetValue( parameters, dataCost, alphasCost, positionCost ) )
    {
    return itk::NumericTraits< double >::max();
    }

  // Add the components
  double  totalCost = dataCost + alphasCost + positionCost;
  if ( std::isnan( totalCost ) || std::isinf( totalCost ) )
    {
    totalCost = itk::NumericTraits< double >::max();
    }

  // Return the result
  return totalCost;

}



//
//
//
bool
AtlasMeshCollectionFastReferencePositionCost
::GetValue( const ParametersType& parameters, double& dataCost, double& alphasCost, double& positionCost ) const
{
  // std::cout << "FastReferencePositionCost: Trying position " << parameters << std::endl;

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
  const double x0 = parameters[ 0 ]; 
  const double y0 = parameters[ 1 ]; 
  const double z0 = parameters[ 2 ];
  for ( std::vector< VertexNeighboringTetrahedronInfo >::const_iterator tetrahedronInfoIt = m_VertexNeighborhood.begin();
         tetrahedronInfoIt != m_VertexNeighborhood.end();
         ++tetrahedronInfoIt )
    {
    const double x1 = ( *tetrahedronInfoIt ).m_X1;
    const double y1 = ( *tetrahedronInfoIt ).m_Y1;
    const double z1 = ( *tetrahedronInfoIt ).m_Z1;
    const double x2 = ( *tetrahedronInfoIt ).m_X2;
    const double y2 = ( *tetrahedronInfoIt ).m_Y2;
    const double z2 = ( *tetrahedronInfoIt ).m_Z2;
    const double x3 = ( *tetrahedronInfoIt ).m_X3;
    const double y3 = ( *tetrahedronInfoIt ).m_Y3;
    const double z3 = ( *tetrahedronInfoIt ).m_Z3;



    // Test if this tethrahedron is left-rotating. If it's not, stop looping over all affected tetrahedra.
    // Do this by calculating the volume of the tetrahedron; this should be positive.
    // In what follows, the matrix Lambda is the Jacobian of the transform from a standarized tetrahedron
    // ( ( 0 0 0 )^T, ( 1 0 0 )^T, ( 0 1 0 )^T, ( 0 0 1 )^T ), which has volume 1/6, to the actual tetrahedron
    const double  lambda11 = -x0 + x1;
    const double  lambda21 = -y0 + y1;
    const double  lambda31 = -z0 + z1;
    const double  lambda12 = -x0 + x2;
    const double  lambda22 = -y0 + y2;
    const double  lambda32 = -z0 + z2;
    const double  lambda13 = -x0 + x3;
    const double  lambda23 = -y0 + y3;
    const double  lambda33 = -z0 + z3;
    const double  volume = ( lambda11 * ( lambda22*lambda33 - lambda32*lambda23 )
                            - lambda12 * ( lambda21*lambda33 - lambda31*lambda23 )
                            + lambda13 * ( lambda21*lambda32 - lambda31*lambda22 ) ) / 6;
    if ( volume <= 0  )
      {
      // std::cout << "FastReferencePositionCost: referencePosition would cause a tetrahedron flip; rejecting this possibility" << std::endl;
      return false;
      }
      
      
    // Also test badness of tetrahedron  
    AtlasMesh::PointType   point0;
    point0[ 0 ] = x0;
    point0[ 1 ] = y0;
    point0[ 2 ] = z0;

    AtlasMesh::PointType   point1;
    point1[ 0 ] = x1;
    point1[ 1 ] = y1;
    point1[ 2 ] = z1;

    AtlasMesh::PointType   point2;
    point2[ 0 ] = x2;
    point2[ 1 ] = y2;
    point2[ 2 ] = z2;

    AtlasMesh::PointType   point3;
    point3[ 0 ] = x3;
    point3[ 1 ] = y3;
    point3[ 2 ] = z3;

    const double badness = TetrahedronRadiusRatio( point0, point1, point2, point3 );

    //// std::cout << "FastReferencePositionCost: badness of tetrahedron "
    //          << tetrahedronInfoIt->m_TetrahedronId << " is " << badness << std::endl;
    // // std::cout << "                      point0: " << point0 << std::endl;
    // // std::cout << "                      point1: " << point1 << std::endl;
    // // std::cout << "                      point2: " << point2 << std::endl;
    // // std::cout << "                      point3: " << point3 << std::endl;
    if ( badness > 10 )
      {
      // std::cout << "FastReferencePositionCost: referencePosition is too bad ("
      //           << badness << "); rejecting this possibility" << std::endl;
      return false;
      }



    } // End loop over all tetrahedra



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
  m_CostCalculationMeshCollection = AtlasMeshCollection::New();
  m_CostCalculationMeshCollection->SetPositions( positions );
  m_CostCalculationMeshCollection->SetPointParameters( pointParameters );
  m_CostCalculationMeshCollection->SetCells( m_Cells );
  m_CostCalculationMeshCollection->SetReferencePosition( m_ReferencePosition );
  m_CostCalculationMeshCollection->SetK( m_K );
  
  
  // Estimate the mesh collection alphas and positions
  m_Estimator->SetInitialMeshCollection( m_CostCalculationMeshCollection );
  m_Estimator->Estimate();
  //meshCollection->Write( "tmp.txt" ); 

 
  // Compute the total model cost
  m_DataAndAlphasCostCalculator->SetMeshCollection( m_CostCalculationMeshCollection );
  m_DataAndAlphasCostCalculator->GetDataCostAndAlphasCost( dataCost, alphasCost );

  positionCost = 0.0;

  return true;
}




} // End namespace kvl
