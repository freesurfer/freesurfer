#include "kvlAtlasMeshBuilder.h"

#include "kvlAtlasMeshCollectionModelLikelihoodCalculator.h"
#include "kvlParameterOrderPowellOptimizer.h"
#include "kvlAtlasMeshCollectionFastReferencePositionCost.h"





namespace kvl
{


//itk::SimpleFastMutexLock  mutex;
AtlasMeshBuilderMutexLock  mutex;


//
//
//
AtlasMeshBuilder
::AtlasMeshBuilder()
{
  m_InitialSize.Fill( 10 );
  m_InitialStiffnesses = std::vector< double >( 1, 0.1 );
  m_Mesher = MultiResolutionAtlasMesher::New();

  m_PowellAbsolutePrecision = 1.0;

  m_StuckCount = 0;

  m_Current = 0;
  m_Progress = 0;
  m_Verbose = false;
  
  m_IterationNumber = 0;
  m_MaximumNumberOfIterations = 30;

  m_RetainedMiniCollection = 0;
  m_CollapsedMiniCollection = 0;

  m_RetainedCost = 0;
  m_RetainedDataCost = 0;
  m_RetainedAlphasCost = 0;
  m_RetainedPositionCost = 0;

  m_CollapsedCost = 0;
  m_CollapsedDataCost = 0;
  m_CollapsedAlphasCost = 0;
  m_CollapsedPositionCost = 0;

  m_EdgeCollapseEncouragementFactor = 1.0;
}



//
//
//
AtlasMeshBuilder
::~AtlasMeshBuilder()
{
}




//
//
//
void
AtlasMeshBuilder
::PrintSelf( std::ostream& os, itk::Indent indent ) const
{



}



//
//
//
void
AtlasMeshBuilder
::SetUp( const std::vector< LabelImageType::ConstPointer >& labelImages,
         const CompressionLookupTable*  compressionLookupTable,
         const itk::Size< 3>&  initialSize, 
         const std::vector< double >& initialStiffnesses,
         const unsigned int maximumNumberOfIterations)
{
  m_LabelImages = labelImages;
  m_CompressionLookupTable = compressionLookupTable;
  m_InitialSize = initialSize;
  m_InitialStiffnesses = initialStiffnesses;
  m_Mesher->SetUp( m_LabelImages, m_CompressionLookupTable, m_InitialSize, m_InitialStiffnesses );
  m_MaximumNumberOfIterations = maximumNumberOfIterations;
}




//
//
//
void
AtlasMeshBuilder
::Build( AtlasMeshCollection* explicitStartCollection, double edgeCollapseEncouragementFactor )
{

  //
  m_EdgeCollapseEncouragementFactor = edgeCollapseEncouragementFactor;
  
  // Estimate high-resolution mesh first
  if ( explicitStartCollection == 0 )
    {
    m_Mesher->Go();
    m_Current = const_cast< AtlasMeshCollection* >( m_Mesher->GetCurrentMeshCollection() );
    }
  else
    {
    if ( fabs( explicitStartCollection->GetK() - m_InitialStiffnesses.back() ) < 1e-2 )
      {
      m_Current = explicitStartCollection;
      }
    else
      {
      // 
      std::cout << "Got an explicitStartCollection, but it's stiffness does not match: " 
                << explicitStartCollection->GetK() << " vs. " << m_InitialStiffnesses.back() << std::endl;
      std::cout << "    difference is " << explicitStartCollection->GetK() - m_InitialStiffnesses.back() << std::endl;          
      std::cout << "So re-estimating that collection first" << std::endl;
      explicitStartCollection->SetK( m_InitialStiffnesses.back() );
      AtlasParameterEstimator::Pointer  estimator = AtlasParameterEstimator::New();
      estimator->SetLabelImages( m_LabelImages, m_CompressionLookupTable );
      estimator->SetInitialMeshCollection( explicitStartCollection );
      estimator->SetPositionOptimizer( AtlasParameterEstimator::LBFGS );    
      estimator->Estimate( true );
      m_Current = const_cast< AtlasMeshCollection* >( estimator->GetCurrentMeshCollection() );
      m_Current->Write( "debug_explicit_after_estimating.txt" );
      }
    }


  // Now start simplifying the result
  this->InvokeEvent( itk::StartEvent() );
  itk::MultiThreader::Pointer  threader = itk::MultiThreader::New();
  int  numberOfThreads = threader->GetNumberOfThreads();
  std::cout << "Using " << numberOfThreads << " threads" << std::endl;

  // ICM-like algorithm
  for ( m_IterationNumber = 0; m_IterationNumber < m_MaximumNumberOfIterations; m_IterationNumber++ )
    {
    std::cout << "Just started one iteration" << std::endl;
    this->DebugOn();

    this->InvokeEvent( itk::IterationEvent() );

    std::cout << "Iteration " << m_IterationNumber << std::endl;

    // Make a container with all the edges to analyze. The threads will remove the first one that's
    // independent of the edges other threads are currently working on, and analyze that. In order to
    // keep track of what edges are or are not independent, we keep score, for every node, if how many
    // threads are currently working on a related edge.
    std::set< AtlasMesh::CellIdentifier >  edges = this->GetRandomizedEdgesAsSet();
    std::map< AtlasMesh::PointIdentifier, int >  pointOccupancies;
    for ( AtlasMesh::PointsContainer::ConstIterator  it = m_Current->GetReferencePosition()->Begin();
          it != m_Current->GetReferencePosition()->End(); ++it )
      {
      pointOccupancies[ it.Index() ] = 0;
      }


    // Do the job
    m_StuckCount = 0; // Reset from previous iterations
    LoadBalancedThreadStruct  str;
    str.m_Builder = this;
    str.m_Edges = edges;
    str.m_PointOccupancies = pointOccupancies;
    threader->SetSingleMethod( this->LoadBalancedThreaderCallback, &str );
    threader->SingleMethodExecute();


    //
    // std::ostringstream  outFileNameStream;
    // outFileNameStream << "builderMeshCollection_ThreadedIteration" << m_IterationNumber << ".txt";
    // m_Current->Write( outFileNameStream.str().c_str() );

    std::cout << "Just finished one iteration" << std::endl;
    //this->DebugOn();
    } // End ICM iterations

  this->InvokeEvent( itk::EndEvent() );


}




//
//
//
void
AtlasMeshBuilder
::GetCurrentDataAndAlphasCost( double& currentDataCost, double& currentAlphasCost ) const
{

  if ( !m_Current )
  {
    currentDataCost = 0.0;
    currentAlphasCost = 0.0;

    return;
  }

  this->GetDataCostAndAlphasCost( m_Current, currentDataCost, currentAlphasCost );

}



//
//
//
double
AtlasMeshBuilder
::GetCurrentPositionCost() const
{
  if ( !m_Current )
  {
    return 0.0;
  }

  return this->GetPositionCost( m_Current );
}




//
//
//
double
AtlasMeshBuilder
::GetCurrentCost() const
{

  if ( !m_Current )
  {
    return 0.0;
  }

  double  dataCost;
  double  alphasCost;
  this->GetDataCostAndAlphasCost( m_Current, dataCost, alphasCost );
  double  positionCost = this->GetPositionCost( m_Current );

  return dataCost + alphasCost + positionCost;
}



//
// 
//
AtlasMeshCollection::Pointer
AtlasMeshBuilder
::TryToCollapseFast( const AtlasMeshCollection* miniCollection,
                     AtlasMesh::CellIdentifier  edgeId,
                     double& miniDataCost, double& miniAlphasCost, double& miniPositionCost,
                     std::set< AtlasMesh::CellIdentifier >&  disappearingCells ) const
{
  //std::cout << "\n\n\n\n\n==================================================" << std::endl;
  //std::cout << "Trying to collapse edge " << edgeId << std::endl;


  // Get the mesh with the edge collapsed.
  AtlasMesh::CellIdentifier  unifiedVertexId;
  AtlasMeshCollection::Pointer  child;
  if ( !miniCollection->GetCollapsed( edgeId, child, disappearingCells, unifiedVertexId ) )
    {
    return nullptr;
    }

  const AtlasMesh::PointIdentifier  unifiedPointId =
             *( child->GetCells()->ElementAt( unifiedVertexId )->PointIdsBegin() );


  // Set up the cost calculator
  kvl::AtlasMeshCollectionFastReferencePositionCost::ParametersType  initialPosition( 3 );
  initialPosition[ 0 ] = child->GetReferencePosition()->ElementAt( unifiedPointId )[ 0 ];
  initialPosition[ 1 ] = child->GetReferencePosition()->ElementAt( unifiedPointId )[ 1 ];
  initialPosition[ 2 ] = child->GetReferencePosition()->ElementAt( unifiedPointId )[ 2 ];



  // Decide whether or not we're going to optimize the reference position --
  // let's not do it if we're inside a uniform image area as it won't improve
  // things while taking valuable time
  bool  optimizeUnifiedPointReferencePosition = false;

  // Look up the label with the highest alpha in the first point of the miniCollection
  AtlasMesh::PointDataContainer::ConstIterator  pointParamIt = miniCollection->GetPointParameters()->Begin();
  int  maximumAlphaLabelNumber = 0;
  double  maximumAlpha = itk::NumericTraits< double >::min();
  for ( unsigned int classNumber = 0; classNumber < pointParamIt.Value().m_Alphas.Size(); classNumber++ )
    {
    if ( pointParamIt.Value().m_Alphas[ classNumber ] > maximumAlpha )
      {
      maximumAlpha = pointParamIt.Value().m_Alphas[ classNumber ];
      maximumAlphaLabelNumber = classNumber;
      }
    }

  // Look at the alphas in each of the points of miniCollection
  const double threshold = 0.90;
  for ( ; pointParamIt != miniCollection->GetPointParameters()->End(); ++pointParamIt )
    {
    if ( pointParamIt.Value().m_Alphas[ maximumAlphaLabelNumber ] < threshold )
      {
      optimizeUnifiedPointReferencePosition = true;
      }
    }

  // Set up cost function
  kvl::AtlasMeshCollectionFastReferencePositionCost::Pointer  costFunction =
                         kvl::AtlasMeshCollectionFastReferencePositionCost::New();
  costFunction->SetLabelImages( m_LabelImages, m_CompressionLookupTable );
  costFunction->SetPointId( unifiedPointId );
  costFunction->SetInitialMeshCollection( child );
  kvl::AtlasMeshCollectionFastReferencePositionCost::ParametersType  optimalPosition = initialPosition;


  if ( !optimizeUnifiedPointReferencePosition )
    {
    if ( m_Verbose )
      {
      std::cout << "NOT optimizing unified reference position" << std::endl;
      }
    }
  else
    {

    // Optimize the position of the unified vertex in the reference position
    const bool canMoveX = child->GetPointParameters()->ElementAt( unifiedPointId ).m_CanMoveX;
    const bool canMoveY = child->GetPointParameters()->ElementAt( unifiedPointId ).m_CanMoveY;
    const bool canMoveZ = child->GetPointParameters()->ElementAt( unifiedPointId ).m_CanMoveZ;
    if ( canMoveX || canMoveY || canMoveZ )
      {
      // Decide what to optimize
      ParameterOrderPowellOptimizer::ParameterOrderType  parameterOrder( 3 );
      int  cumulativeSum = 0;
      if ( canMoveX )
        {
        parameterOrder[ 0 ] = 1;
        cumulativeSum++;
        }
      else
        {
        parameterOrder[ 0 ] = 0;
        }

      if ( canMoveY )
        {
        parameterOrder[ 1 ] = cumulativeSum + 1;
        cumulativeSum++;
        }
      else
        {
        parameterOrder[ 1 ] = 0;
        }

      if ( canMoveZ )
        {
        parameterOrder[ 2 ] = cumulativeSum + 1;
        }
      else
        {
        parameterOrder[ 2 ] = 0;
        }

      //std::cout << "parameterOrder: " << parameterOrder << std::endl;


      // Optimize the reference position
      //std::cout << "Optimizing the position of the unified point " << unifiedPointId << std::endl;
      ParameterOrderPowellOptimizer::Pointer  optimizer = ParameterOrderPowellOptimizer::New();
      optimizer->SetCostFunction( costFunction );
      optimizer->SetInitialPosition( initialPosition );
      // optimizer->SetAbsolutePrecisionBrent( m_PowellAbsolutePrecision );
      optimizer->SetStepTolerance( m_PowellAbsolutePrecision );
      // optimizer->SetAbsolutePrecision( m_PowellAbsolutePrecision );
      optimizer->SetParameterOrder( parameterOrder );
      optimizer->StartOptimization();

      // Retrieve the optimal reference position
      if ( optimizer->GetCurrentCost() != itk::NumericTraits< double >::max() )
        {
        optimalPosition = optimizer->GetCurrentPosition();
        //std::cout << "Changed reference position for the unified point " << unifiedPointId
        //          << " from " << initialPosition << " to " << optimalPosition << std::endl;
        }
      else
        {
        return nullptr;
        }

      }  // End test if point can move

    } // End test if we are going to optimize position


  // Get the cost at the optimal position
  //std::cout << "Getting cost at optimal position " << optimalPosition << ": ";
  if ( !costFunction->GetValue( optimalPosition, miniDataCost, miniAlphasCost, miniPositionCost ) )
    {
    //std::cout << "not possible" << std::endl;
    return nullptr;
    }
  //std::cout << "       miniDataCost     : " << miniDataCost << std::endl;
  //std::cout << "       miniAlphasCost   : " << miniAlphasCost << std::endl;
  //std::cout << "       miniPositionCost : " << miniPositionCost << std::endl;


  // Return
  if ( std::isnan( miniDataCost + miniAlphasCost + miniPositionCost ) ||
       std::isinf( miniDataCost + miniAlphasCost + miniPositionCost ) )
    {
    return nullptr;
    }
  else
    {
    return const_cast< AtlasMeshCollection* >( costFunction->GetCostCalculationMeshCollection() );
    }

}




//
//
//
AtlasMeshCollection::Pointer
AtlasMeshBuilder
::TryToRetainFast( const AtlasMeshCollection* miniCollectionConst,
                   AtlasMesh::CellIdentifier  edgeId,
                   double& miniDataCost, double& miniAlphasCost, double& miniPositionCost )
{

  //std::cout << "\n\n\n\n\n==================================================" << std::endl;
  //std::cout << "Trying to retain edge " << edgeId << std::endl;

  AtlasMeshCollection::Pointer  miniCollection = this->GetFakeCopy( miniCollectionConst );


  // Retrieve the point ids of the vertices of this edge
  AtlasMesh::CellType::PointIdConstIterator  pointIt =
                 miniCollection->GetCells()->ElementAt( edgeId )->PointIdsBegin();
  const AtlasMesh::PointIdentifier  point0Id = *pointIt;
  ++pointIt;
  const AtlasMesh::PointIdentifier  point1Id = *pointIt;


  // Retrieve the vertex ids of the vertices of this edge, by looping over all vertex cells and
  // checking their point id
  AtlasMesh::CellIdentifier  vertex0Id;
  AtlasMesh::CellIdentifier  vertex1Id;
  for ( AtlasMesh::CellsContainer::ConstIterator  cellIt = miniCollection->GetCells()->Begin();
        cellIt != miniCollection->GetCells()->End(); ++cellIt )
    {
    const AtlasMesh::CellType*  cell = cellIt.Value();
    if( cell->GetType() != AtlasMesh::CellType::VERTEX_CELL )
      {
      continue;
      }

    if ( *( cell->PointIdsBegin() ) == point0Id )
      {
      vertex0Id = cellIt.Index();
      //std::cout << "vertex0Id: " << vertex0Id << std::endl;
      }
    else if (  *( cell->PointIdsBegin() ) == point1Id )
      {
      vertex1Id = cellIt.Index();
      //std::cout << "vertex1Id: " << vertex1Id << std::endl;
      }

    }


  // Optimize position of vertices.
  // TODO: should be iterated!
  if ( !this->OptimizeReferencePositionFast( miniCollection, vertex0Id, false ) )
    {
    return nullptr;
    }
  if ( !this->OptimizeReferencePositionFast( miniCollection, vertex1Id, false ) )
    {
    return nullptr;
    }


  // Measure the cost after optimization
  this->GetDataCostAndAlphasCost( miniCollection, miniDataCost, miniAlphasCost );
  miniPositionCost = this->GetPositionCost( miniCollection );


  // Return
  if ( ( std::isnan( miniDataCost + miniAlphasCost + miniPositionCost ) ) ||
       ( std::isinf( miniDataCost + miniAlphasCost + miniPositionCost ) ) )
    {
    return nullptr;
    }
  else
    {
    return miniCollection;
    }


}



//
//
//
bool
AtlasMeshBuilder
::OptimizeReferencePositionFast( AtlasMeshCollection* meshCollection,
                                 AtlasMesh::CellIdentifier  vertexId, bool optimize ) const
{

  // Retrieve the point id of the point corresponding to the vertex
  const AtlasMesh::PointIdentifier  pointId =
          *( meshCollection->GetCells()->ElementAt( vertexId )->PointIdsBegin() );

  // Set up the cost calculator
  kvl::AtlasMeshCollectionFastReferencePositionCost::ParametersType  initialPosition( 3 );
  initialPosition[ 0 ] = meshCollection->GetReferencePosition()->ElementAt( pointId )[ 0 ];
  initialPosition[ 1 ] = meshCollection->GetReferencePosition()->ElementAt( pointId )[ 1 ];
  initialPosition[ 2 ] = meshCollection->GetReferencePosition()->ElementAt( pointId )[ 2 ];
  kvl::AtlasMeshCollectionFastReferencePositionCost::Pointer  costFunction =
                                kvl::AtlasMeshCollectionFastReferencePositionCost::New();
  costFunction->SetLabelImages( m_LabelImages, m_CompressionLookupTable );
  costFunction->SetPointId( pointId );
  costFunction->SetInitialMeshCollection( meshCollection );
  kvl::AtlasMeshCollectionFastReferencePositionCost::ParametersType  optimalPosition = initialPosition;


  // Optimize the position of the vertex in the reference position
  const bool canMoveX = meshCollection->GetPointParameters()->ElementAt( pointId ).m_CanMoveX;
  const bool canMoveY = meshCollection->GetPointParameters()->ElementAt( pointId ).m_CanMoveY;
  const bool canMoveZ = meshCollection->GetPointParameters()->ElementAt( pointId ).m_CanMoveZ;
  if ( !canMoveX && !canMoveY && !canMoveZ )
    {
    return true;
    }

  if ( optimize )
    {
    // Decide what to optimize
    ParameterOrderPowellOptimizer::ParameterOrderType  parameterOrder( 3 );
    int  cumulativeSum = 0;
    if ( canMoveX )
      {
      parameterOrder[ 0 ] = 1;
      cumulativeSum++;
      }
    else
      {
      parameterOrder[ 0 ] = 0;
      }

    if ( canMoveY )
      {
      parameterOrder[ 1 ] = cumulativeSum + 1;
      cumulativeSum++;
      }
    else
      {
      parameterOrder[ 1 ] = 0;
      }

    if ( canMoveZ )
      {
      parameterOrder[ 2 ] = cumulativeSum + 1;
      }
    else
      {
      parameterOrder[ 2 ] = 0;
      }

    // std::cout << "parameterOrder: " << parameterOrder << std::endl;


    // Optimize the reference position
    // std::cout << "Optimizing the position of the unified point " << pointId << std::endl;
    ParameterOrderPowellOptimizer::Pointer  optimizer = ParameterOrderPowellOptimizer::New();
    optimizer->SetCostFunction( costFunction );
    optimizer->SetInitialPosition( initialPosition );
    // optimizer->SetAbsolutePrecisionBrent( m_PowellAbsolutePrecision );
    optimizer->SetStepTolerance( m_PowellAbsolutePrecision );
    // optimizer->SetAbsolutePrecision( m_PowellAbsolutePrecision );
    optimizer->SetParameterOrder( parameterOrder );
    optimizer->StartOptimization();

    // Retrieve the optimal reference position
    if ( optimizer->GetCurrentCost() != itk::NumericTraits< double >::max() )
      {
      optimalPosition = optimizer->GetCurrentPosition();
      // std::cout << "                               Changed reference position for the point " << pointId
      //           << " from " << initialPosition << " to " << optimalPosition << std::endl;
      // const double  initialCost = costFunction->GetValue( initialPosition );
      // const double  optimalCost = costFunction->GetValue( optimalPosition );
      // std::cout << "                                   " << initialCost << "  ->  " << optimalCost << std::endl;
      }
    else
      {
      return false;
      }

    } // End test if we need to optimize


  // Get the optimized outer mini mesh from the costFunction
  if ( costFunction->GetValue( optimalPosition ) == itk::NumericTraits< double >::max() )
    {
    return false;
    }
  AtlasMeshCollection::Pointer  optimizedCollection =
    const_cast< AtlasMeshCollection* >( costFunction->GetCostCalculationMeshCollection() );

  meshCollection->SetPositions( optimizedCollection->GetPositions() );
  meshCollection->SetReferencePosition( optimizedCollection->GetReferencePosition() );
  meshCollection->SetK( meshCollection->GetK() );
  meshCollection->SetPointParameters( optimizedCollection->GetPointParameters() );

  return true;
}





//
//
//
double
AtlasMeshBuilder
::GetPositionCost( const AtlasMeshCollection* meshCollection ) const
{
  return 0.0;
}




//
//
//
void
AtlasMeshBuilder
::GetDataCostAndAlphasCost( const AtlasMeshCollection* meshCollection, double& dataCost, double& alphasCost ) const
{

  AtlasMeshCollectionModelLikelihoodCalculator::Pointer  calculator = 
                                         AtlasMeshCollectionModelLikelihoodCalculator::New();
  calculator->SetMeshCollection( meshCollection );
  calculator->SetLabelImages( m_LabelImages, m_CompressionLookupTable );

  calculator->GetDataCostAndAlphasCost( dataCost, alphasCost );

}



//
//
//
std::vector< AtlasMesh::CellIdentifier >
AtlasMeshBuilder
::Permutate(  const std::vector< AtlasMesh::CellIdentifier >&  edgeList ) const
{

  // Populate sortable vector
  std::vector< EdgeElement >    edgeElements;
  for ( std::vector< AtlasMesh::CellIdentifier >::const_iterator  it = edgeList.begin(); it != edgeList.end(); ++it )
    {
    edgeElements.push_back( EdgeElement( *it ) );
    }

  // Sort the sortable vector
  std::sort( edgeElements.begin(), edgeElements.end() );

  // Construct something to return
  std::vector< AtlasMesh::CellIdentifier >  result;
  for ( std::vector< EdgeElement >::const_iterator  it = edgeElements.begin(); it != edgeElements.end(); ++it )
    {
    result.push_back( ( *it ).GetEdgeId()  );
    }
  return result;
}



//
//
//
AtlasMeshCollection::Pointer
AtlasMeshBuilder
::GetFakeCopy( const AtlasMeshCollection*  input ) const
{
  // Copy the reference position
  AtlasMesh::PointsContainer::Pointer  outputReferencePosition = AtlasMesh::PointsContainer::New();
  for ( AtlasMesh::PointsContainer::ConstIterator  it = input->GetReferencePosition()->Begin();
        it != input->GetReferencePosition()->End(); ++it )
    {
    outputReferencePosition->InsertElement( it.Index(), it.Value() );
    }

  // Copy the points data
  AtlasMesh::PointDataContainer::Pointer  outputPointParameters = AtlasMesh::PointDataContainer::New();
  for ( AtlasMesh::PointDataContainer::ConstIterator  it = input->GetPointParameters()->Begin();
        it != input->GetPointParameters()->End(); ++it )
    {
    outputPointParameters->InsertElement( it.Index(), it.Value() );
    }

  // Copy the positions
  std::vector< AtlasMesh::PointsContainer::Pointer >  outputPositions;
  for ( unsigned int meshNumber = 0; meshNumber < input->GetNumberOfMeshes(); meshNumber++ )
    {
    AtlasMesh::PointsContainer::ConstPointer  inputPosition = ( input->GetPositions()[ meshNumber ] ).GetPointer();
    AtlasMesh::PointsContainer::Pointer  outputPosition = AtlasMesh::PointsContainer::New();
    for ( AtlasMesh::PointsContainer::ConstIterator  it = inputPosition->Begin(); it != inputPosition->End(); ++it  )
      {
      outputPosition->InsertElement( it.Index(), it.Value() );
      }
    outputPositions.push_back( outputPosition );
    }


  // Collect everything and return
  AtlasMeshCollection::Pointer  output = AtlasMeshCollection::New();
  output->SetPointParameters( outputPointParameters );
  output->SetCells( const_cast< AtlasMesh::CellsContainer* >( input->GetCells() ) );
  output->SetReferencePosition( outputReferencePosition );
  output->SetPositions( outputPositions );
  output->SetK( input->GetK() );

  return output;
}





//
//
//
bool
AtlasMeshBuilder
::LoadBalancedAnalyzeEdgeFast( std::set< AtlasMesh::CellIdentifier >&  edges,
                                std::map< AtlasMesh::PointIdentifier, int >&  pointOccupancies,
                                int  threadId )
{

#if 0  
  mutex.Lock(); 
#else
  AtlasMeshBuilderHelper  helper( mutex, pointOccupancies );
#endif
  //std::ostringstream  descriptionStream;
  //descriptionStream << "    [THREAD " << threadId << "] locked mutex because trying to find edge to analyze";
  //mutex.DescriptiveLock( descriptionStream.str() );

  // Check if there is anything left to do
  if ( edges.size() == 0 )
    {
#if 0      
    mutex.Unlock();
    //mutex.DescriptiveUnlock();
#endif    
    return false;
    }


  // Look for the first edge that can be safely worked on
  // std::cout << "    [THREAD " << threadId << "] " << " looking for the first edge that can be safely worked on" << std::endl;
  AtlasMesh::CellIdentifier  edgeId = itk::NumericTraits< AtlasMesh::CellIdentifier >::max();
  std::vector< AtlasMesh::CellIdentifier >  nonExistingEdgesEncountered;
  for ( std::set< AtlasMesh::CellIdentifier >::const_iterator  it = edges.begin();
        it != edges.end(); ++it )
    {
    // Test if the edge actually still exists in the mesh. If not, let's tag it for removal
    if ( !m_Current->GetCells()->IndexExists( *it ) )
      {
      nonExistingEdgesEncountered.push_back( *it );
      continue;
      }


    // Test if this edge can be safely worked on
    AtlasMesh::CellType::PointIdConstIterator  pit = m_Current->GetCells()->ElementAt( *it )->PointIdsBegin();
    const AtlasMesh::PointIdentifier  p0Id = *pit;
    ++pit;
    const AtlasMesh::PointIdentifier  p1Id = *pit;

    if ( ( pointOccupancies[ p0Id ] == 0 ) && ( pointOccupancies[ p1Id ] == 0 ) )
      {
      // OK, this edge is safe. Remember it's id.
      edgeId = *it;
      break;
      }

    //std::cout << "    [THREAD " << threadId << "] " << " Hmm... edge " << *it << " is protected" << std::endl;
    } // End loop over all edges


  // Remove the non-existing edges we're encountered
  for ( std::vector< AtlasMesh::CellIdentifier >::const_iterator  it = nonExistingEdgesEncountered.begin();
        it != nonExistingEdgesEncountered.end(); ++it )
    {
    if ( m_Verbose )
      {
      std::cout << "    [THREAD " << threadId << "] " << "Encountered non-existing edge " << *it << ". Erasing it" << std::endl;
      }
    edges.erase( *it );
    }


  if ( edgeId == itk::NumericTraits< AtlasMesh::CellIdentifier >::max() )
    {
    // There were edges to be analyzed, but can't work on them at this time cause others or
    // working on it
    if ( m_Verbose )
      {
      std::cout << "    [THREAD " << threadId << "] " << "There are still " << edges.size()
                << " edges but they are all protected at this point" << std::endl;
      }          

    m_StuckCount++; // This is not actually thread-safe, but we don't care if this count goes up a bit
                    // slower than you'd expect (in the unlikely event that two threads try to increase 
                    // this count at exactly the same time)
    const double  averageNumberOfSecondsToSleep = 60.0;
    const int  stuckCountToGiveUpAt = 30 * 20;  // If every thread sleeps 1min on average, and 20 threads
                                                // are used, this will give up after 30min
    if ( m_StuckCount % 10 == 0 )
      {  
      std::cout << "    m_StuckCount is at " << m_StuckCount 
                << " (exiting at " << stuckCountToGiveUpAt << ")" << std::endl; 
      }

    if ( m_StuckCount >= stuckCountToGiveUpAt )
      {
      std::ostringstream  stuckStream;
      stuckStream << "threadsNeverReturning";
      m_Current->Write( stuckStream.str().c_str() );
      exit(-1);
      //itkExceptionMacro( << "Error: threads never returning" );
      }

#if 0      
    mutex.Unlock();
    //mutex.DescriptiveUnlock();
#else
    helper.Unlock();
#endif
    
    // Sleep averageNumberOfSecondsToSleep (on average)
    const int  maxRange = ( 1000000 * averageNumberOfSecondsToSleep * 2 ); // x2 because of average
    const int  numberOfMicrosecondsToSleep = rand() % maxRange; // Integer between 0 and (maxRange-1)
    usleep( numberOfMicrosecondsToSleep );
    
    return true;
    }


  AtlasMesh::CellType::PointIdConstIterator  pit
                 = m_Current->GetCells()->ElementAt( edgeId )->PointIdsBegin();
  const AtlasMesh::PointIdentifier  p0Id = *pit;
  ++pit;
  const AtlasMesh::PointIdentifier  p1Id = *pit;
  if ( m_Verbose )
    {
    std::cout << "    [THREAD " << threadId << "] " << "    Analyzing edge with id: " << edgeId
              << " (pointIds " << p0Id << " and " << p1Id << " ) "
              << " (reference positions " << m_Current->GetReferencePosition()->ElementAt( p0Id )
              << " and " << m_Current->GetReferencePosition()->ElementAt( p1Id ) << ")" << std::endl;
    }          


  // Remove this edge from the edges to be analyzed
  if ( m_Verbose )
    {
    std::cout << "    [THREAD " << threadId << "] " 
              << "         Removing edge with id: " << edgeId 
              << " from edges (" << edges.size() << ")" << std::endl;
    }
  edges.erase( edgeId );

    
  // Get the mini collection surrounding the edge to try to collapse
  AtlasMeshCollection::ConstPointer  miniCollection = m_Current->GetRegionGrown( edgeId, 1 ).GetPointer();


  // Collect the points that are affected by our edge
  std::vector< AtlasMesh::PointIdentifier >  affectedPoints;
  for ( AtlasMesh::PointsContainer::ConstIterator  refPosIt = miniCollection->GetReferencePosition()->Begin();
        refPosIt != miniCollection->GetReferencePosition()->End(); ++refPosIt )
    {
    // std::set's insert only inserts if element doesn't exist already, so no worries about double entries
    affectedPoints.push_back( refPosIt.Index() );
    }

#if 0    
  // Indicate that we are working on them
  for ( std::vector< AtlasMesh::PointIdentifier >::const_iterator  it = affectedPoints.begin();
        it != affectedPoints.end(); ++it )
    {
    //std::cout << "    [THREAD " << threadId << "] " << "        Edge " << edgeId << " is now locking point " << *it << std::endl;
    ( pointOccupancies[ *it ] )++;
    }

  mutex.Unlock();
  //mutex.DescriptiveUnlock();
#else
  helper.SetAffectedPoints( affectedPoints );
  helper.Unlock();
#endif
  
  // itk::TimeProbe  timeProbe;
  // timeProbe.Start();

  //  Calculate cost when you just optimize this edge's vertex positions
  if ( m_Verbose )
    {
    std::cout << "    [THREAD " << threadId << "] " << "Trying to retain edge " << edgeId << std::endl;
    }
  double  retainedCost = 0;
  double  retainedDataCost = 0;
  double  retainedAlphasCost = 0;
  double  retainedPositionCost = 0;
  AtlasMeshCollection::Pointer  retainedMiniCollection =
     this->TryToRetainFast( miniCollection, edgeId,
                            retainedDataCost, retainedAlphasCost, retainedPositionCost );
  if ( !retainedMiniCollection )
    {
    // Couldn't retain this edge
    retainedCost = itk::NumericTraits< double >::max();
    }
  else
    {
    retainedCost = retainedDataCost + retainedAlphasCost + retainedPositionCost;
    }



  // Calculate the cost of an edge collapse
  if ( m_Verbose )
    {
    std::cout << "    [THREAD " << threadId << "] " << "Trying to collapse edge " << edgeId << std::endl;
    }
  double  collapsedCost = 0;
  double  collapsedDataCost = 0;
  double  collapsedAlphasCost = 0;
  double  collapsedPositionCost = 0;
  std::set< AtlasMesh::CellIdentifier >  collapsedDisappearingCells;
  AtlasMeshCollection::Pointer  collapsedMiniCollection =
     this->TryToCollapseFast( miniCollection, edgeId,
                              collapsedDataCost, collapsedAlphasCost, collapsedPositionCost,
                              collapsedDisappearingCells );
  if ( !collapsedMiniCollection )
    {
    // Couldn't collapse this edge.
    collapsedCost = itk::NumericTraits< double >::max();
    }
  else
    {
    collapsedCost = collapsedDataCost + collapsedAlphasCost + collapsedPositionCost;
    }


  // Evaluate which move is best
  if ( m_Verbose )
    {
    std::cout << "    [THREAD " << threadId << "] " << "                 retainedCost : " << retainedCost
              << "  (" << retainedDataCost << " + " << retainedAlphasCost
              << " + " << retainedPositionCost <<") " << std::endl;
    std::cout << "    [THREAD " << threadId << "] " << "                 collapsedCost: " << collapsedCost
              << "  (" << collapsedDataCost << " + " << collapsedAlphasCost
              << " + " << collapsedPositionCost <<") " << std::endl;
    }          
  std::vector< double >  totalCosts;
  totalCosts.push_back( retainedCost * m_EdgeCollapseEncouragementFactor );
  totalCosts.push_back( collapsedCost );
  double  minTotalCost = itk::NumericTraits< double >::max();
  int minTotalCostIndex = -1;
  for ( unsigned int i = 0; i < totalCosts.size(); i++ )
    {
    if ( totalCosts[ i ] < minTotalCost )
      {
      minTotalCost = totalCosts[ i ];
      minTotalCostIndex = i;
      }
    }


  // timeProbe.Stop();
  // std::cout << "    [THREAD " << threadId << "] took " << timeProbe.GetMeanTime() << " seconds to evaluate moves" << std::endl;


  if ( minTotalCostIndex == -1 )
    {
#if 0      
    mutex.Lock();
#else
    helper.Lock();
#endif
    // descriptionStream << "    [THREAD " << threadId << "] locked mutex because impossible configuration encountered";
    // std::ostringstream  descriptionStream;
    // mutex.DescriptiveLock( descriptionStream.str() );

    std::cout << "    [THREAD " << threadId << "] " << "Impossible configuration encountered at eget with id " << edgeId << std::endl;
    //std::ostringstream  impossibleStream;
    //impossibleStream << "impossible_" << edgeId;
    //m_Current->Write( impossibleStream.str().c_str() );
    //minTotalCostIndex = 0;

#if 0    
    // Now "un-protect" the points we flagged as being worked on
    for ( std::vector< AtlasMesh::PointIdentifier >::const_iterator  it = affectedPoints.begin();
          it != affectedPoints.end(); ++it )
      {
      //std::cout << "    [THREAD " << threadId << "] " << "        Edge " << edgeId << " is now unlocking point " << *it << std::endl;
      ( pointOccupancies[ *it ] )--;
      }

    mutex.Unlock();
    //mutex.DescriptiveUnlock();
#endif
    
    return true;
    }


#if 0    
  mutex.Lock();
  // std::ostringstream  descriptionStream;
  // descriptionStream << "    [THREAD " << threadId << "] locked mutex because applying best move to the global mesh";
  // mutex.DescriptiveLock( descriptionStream.str() );
#else
  helper.Lock();
#endif
  
  // Do the best move
  if ( minTotalCostIndex == 0 )
    {
    if ( m_Verbose )
      {
      std::cout << "    [THREAD " << threadId << "] " << "        => retaining edge is best solution" << std::endl;
      }

    // Look up the ids of the two points on the edge
    AtlasMesh::CellType::PointIdConstIterator  pointIt = m_Current->GetCells()->ElementAt( edgeId )->PointIdsBegin();
    const AtlasMesh::PointIdentifier  edgePoint0Id = *pointIt;
    ++pointIt;
    const AtlasMesh::PointIdentifier  edgePoint1Id = *pointIt;

    // Copy the reference position of each of the two points from the min mesh to the full mesh
    m_Current->GetReferencePosition()->ElementAt( edgePoint0Id ) =
      retainedMiniCollection->GetReferencePosition()->ElementAt( edgePoint0Id );
    m_Current->GetReferencePosition()->ElementAt( edgePoint1Id ) =
      retainedMiniCollection->GetReferencePosition()->ElementAt( edgePoint1Id );

    // Copy the positions of each of the two points from the mini mesh to the full mesh
    for ( unsigned int meshNumber = 0; meshNumber < m_Current->GetNumberOfMeshes(); meshNumber++ )
      {
      m_Current->GetPositions()[ meshNumber ]->ElementAt( edgePoint0Id ) =
        retainedMiniCollection->GetPositions()[ meshNumber ]->ElementAt( edgePoint0Id );
      m_Current->GetPositions()[ meshNumber ]->ElementAt( edgePoint1Id ) =
        retainedMiniCollection->GetPositions()[ meshNumber ]->ElementAt( edgePoint1Id );
      }

    // Copy the point parameters of each of the two points from the mini mesh to the full mesh
    m_Current->GetPointParameters()->ElementAt( edgePoint0Id ) =
      retainedMiniCollection->GetPointParameters()->ElementAt( edgePoint0Id );
    m_Current->GetPointParameters()->ElementAt( edgePoint1Id ) =
      retainedMiniCollection->GetPointParameters()->ElementAt( edgePoint1Id );

    // You manually modified the reference position; make sure the ReferenceTetrahedronInfos is modified manually
    // as well so things are up-to-date
    for ( AtlasMesh::CellDataContainer::ConstIterator infoIt =
            retainedMiniCollection->GetReferenceTetrahedronInfos()->Begin();
          infoIt != retainedMiniCollection->GetReferenceTetrahedronInfos()->End(); ++infoIt )
      {
      const_cast< AtlasMesh::CellDataContainer* >( m_Current->GetReferenceTetrahedronInfos() )->ElementAt( infoIt.Index() )
      = infoIt.Value();
      }

    }
  else if ( minTotalCostIndex == 1 )
    {     
    if ( m_Verbose )
      {
      std::cout << "    [THREAD " << threadId << "] " << "        => collapsing edge is best solution" << std::endl;
      }
    
    // Look up the ids of the two points on the edge
    AtlasMesh::CellType::PointIdConstIterator  pointIt = m_Current->GetCells()->ElementAt( edgeId )->PointIdsBegin();
    const AtlasMesh::PointIdentifier  edgePoint0Id = *pointIt;
    ++pointIt;
    const AtlasMesh::PointIdentifier  edgePoint1Id = *pointIt;


    // Let's quickly collapse a little bit bigger mesh, and get its cell links. The cell links belonging to the
    // inner points, i.e. the points of our miniCollection, will be copied later on
    AtlasMeshCollection::ConstPointer  biggerMiniCollection = m_Current->GetRegionGrown( edgeId, 2 ).GetPointer();
    if ( !biggerMiniCollection )
      {
      std::cout << "    [THREAD " << threadId << "] " << "Ouch!" << std::endl;
      //exit( -1 );
      itkExceptionMacro( << "Ouch!" );
      }
    AtlasMeshCollection::Pointer  collapsedBiggerMiniCollection;
    std::set< AtlasMesh::CellIdentifier >  dummyDisappearingCells;
    AtlasMesh::CellIdentifier  dummyUnifiedVertexId;
    biggerMiniCollection->GetCollapsed( edgeId, collapsedBiggerMiniCollection,
                                        dummyDisappearingCells, dummyUnifiedVertexId );
    AtlasMesh::CellLinksContainerPointer  collapsedBiggerMiniCollectionCellLinks = collapsedBiggerMiniCollection->GetCellLinks();


    // Do the basic mesh surgery: loop over all cells in the original mini mesh, and look for them
    // in the collapsed mini mesh. If it is found there, delete the corresponding cell in the global mesh,
    // put its pointer to the cell in the collapsed mini mesh, and remove the cell from the collapsed mini mesh
    // (the latter step is necessary as the destructor of the collapsed mini mesh will call delete on the pointer).
    // If it is not found, the cell has disappeared due to the edge collapse, so delete the corresponding cell from
    // the global mesh
    for ( AtlasMesh::CellsContainer::ConstIterator  cellIt = miniCollection->GetCells()->Begin();
          cellIt != miniCollection->GetCells()->End(); ++cellIt )
      {
      if ( collapsedMiniCollection->GetCells()->IndexExists( cellIt.Index() ) )
        {
        // Delete the cell object to which we're currently pointing
        delete m_Current->GetCells()->ElementAt( cellIt.Index() );

        // Instead point to the (corrected) cell object in the collapsed mini mesh
        m_Current->GetCells()->ElementAt( cellIt.Index() ) =
          collapsedMiniCollection->GetCells()->ElementAt( cellIt.Index() );

        // Make sure the destructor of the collapsed mini mesh won't find this cell
        // and delete it!
        collapsedMiniCollection->GetCells()->DeleteIndex( cellIt.Index() );
        }
      else
        {
        // Delete the cell object to which we're currently pointing
        delete m_Current->GetCells()->ElementAt( cellIt.Index() );

        // Remove the index from the container
        m_Current->GetCells()->DeleteIndex( cellIt.Index() );

        // If this cell is a disappearing tetrahedron, also remove it from the ReferenceTetrahedronInfos.
        // We don't have to test if it actually exists or not: std::map::erase() is fool-proof
        const_cast< AtlasMesh::CellDataContainer* >( m_Current->GetReferenceTetrahedronInfos() )->DeleteIndex( cellIt.Index() );
        }

      } // End loop over all cells that need removal or surgery


    // We've altered the mesh connnectivity, so make sure the cell links are up-to-date by copying them
    // from the collapsed mesh. Note: the second point on the edge is actually disappearing; we'll take
    // take care of that later
    for ( AtlasMesh::PointsContainer::ConstIterator pointIt = miniCollection->GetReferencePosition()->Begin();
          pointIt != miniCollection->GetReferencePosition()->End(); ++ pointIt )
      {
      m_Current->GetCellLinks()->ElementAt( pointIt.Index() ) =
        collapsedBiggerMiniCollectionCellLinks->ElementAt( pointIt.Index() );

      // std::cout << "    [THREAD " << threadId << "] " << pointIt.Index() << ": copied "
      //           << m_Current->GetCellLinks()->ElementAt( pointIt.Index() ).size()
      //           << " cell indices from mini " << std::endl;
      }



    // Delete the second point (edgePoint1Id) from the point containers
    m_Current->GetReferencePosition()->DeleteIndex( edgePoint1Id );
    for ( unsigned int meshNumber = 0; meshNumber < m_Current->GetNumberOfMeshes(); meshNumber++ )
      {
      m_Current->GetPositions()[ meshNumber ]->DeleteIndex( edgePoint1Id );
      }
    m_Current->GetPointParameters()->DeleteIndex( edgePoint1Id );
    m_Current->GetCellLinks()->DeleteIndex( edgePoint1Id );


    // Copy the reference position of the unified point from the mini collection to the global collection
    m_Current->GetReferencePosition()->ElementAt( edgePoint0Id ) =
      collapsedMiniCollection->GetReferencePosition()->ElementAt( edgePoint0Id );

    // Copy the positions of the unified point from the mini collection to the global collection
    for ( unsigned int meshNumber = 0; meshNumber < m_Current->GetNumberOfMeshes(); meshNumber++ )
      {
      m_Current->GetPositions()[ meshNumber ]->ElementAt( edgePoint0Id ) =
        collapsedMiniCollection->GetPositions()[ meshNumber ]->ElementAt( edgePoint0Id );
      }

    // Copy the point parameters of the unified point from the mini collection to the global collection
    m_Current->GetPointParameters()->ElementAt( edgePoint0Id ) =
      collapsedMiniCollection->GetPointParameters()->ElementAt( edgePoint0Id );


    // You manually modified the reference position; make sure the ReferenceTetrahedronInfos is modified manually
    // as well so things are up-to-date
    for ( AtlasMesh::CellDataContainer::ConstIterator infoIt =
            collapsedMiniCollection->GetReferenceTetrahedronInfos()->Begin();
          infoIt != collapsedMiniCollection->GetReferenceTetrahedronInfos()->End(); ++infoIt )
      {
      const_cast< AtlasMesh::CellDataContainer* >( m_Current->GetReferenceTetrahedronInfos() )->ElementAt( infoIt.Index() )
      = infoIt.Value();
      }

    } // End deciding what the best operation is for this edge


#if 0    
  // Now "un-protect" the points we flagged as being worked on
  for ( std::vector< AtlasMesh::PointIdentifier >::const_iterator  it = affectedPoints.begin();
        it != affectedPoints.end(); ++it )
    {
    //std::cout << "    [THREAD " << threadId << "] " << "        Edge " << edgeId << " is now unlocking point " << *it << std::endl;
    ( pointOccupancies[ *it ] )--;
    }

  if ( m_Verbose )
    {
    std::cout << "    [THREAD " << threadId << "] " << "Done with edge " << edgeId << std::endl;
    }
  
  mutex.Unlock();
  //mutex.DescriptiveUnlock();
#endif
  
  return true;

}




//
//
//
std::vector< AtlasMesh::CellIdentifier >
AtlasMeshBuilder
::GetRandomizedEdges() const
{

  // Cache edges at this point
  std::vector< AtlasMesh::CellIdentifier >  edges;
  for ( AtlasMesh::CellsContainer::ConstIterator  cellIt = m_Current->GetCells()->Begin();
        cellIt != m_Current->GetCells()->End(); ++cellIt )
    {
    if ( cellIt.Value()->GetType() == AtlasMesh::CellType::LINE_CELL )
      {
      edges.push_back( cellIt.Index() );
      }
    }

  // Randomize the edges' order
  edges =  this->Permutate(  edges );

  return edges;
}



//
//
//
std::set< AtlasMesh::CellIdentifier >
AtlasMeshBuilder
::GetRandomizedEdgesAsSet() const
{

  // Convert from vector to set
  std::vector< AtlasMesh::CellIdentifier >  edgesVector = this->GetRandomizedEdges();
  std::set< AtlasMesh::CellIdentifier >  edges;
  for ( std::vector< AtlasMesh::CellIdentifier >::const_iterator  it = edgesVector.begin();
        it != edgesVector.end(); ++it )
    {
    edges.insert( *it );
    }

  return edges;
}



//
//
//
#if ITK_VERSION_MAJOR >= 5
itk::ITK_THREAD_RETURN_TYPE
AtlasMeshBuilder
::LoadBalancedThreaderCallback( void *arg )
#else  
ITK_THREAD_RETURN_TYPE
AtlasMeshBuilder
::LoadBalancedThreaderCallback( void *arg )
#endif    
{

  // Retrieve the input arguments
#if ITK_VERSION_MAJOR >= 5
  const int  threadId = ((itk::MultiThreaderBase::WorkUnitInfo *)(arg))->WorkUnitID;
  //const int  threadCount = ((itk::MultiThreaderBase::WorkUnitInfo *)(arg))->NumberOfWorkUnits;

  LoadBalancedThreadStruct*  str = (LoadBalancedThreadStruct *)(((itk::MultiThreaderBase::WorkUnitInfo *)(arg))->UserData);  
#else  
  const int  threadId = ((itk::MultiThreader::ThreadInfoStruct *)(arg))->ThreadID;
  //const int  threadCount = ((itk::MultiThreader::ThreadInfoStruct *)(arg))->NumberOfThreads;

  LoadBalancedThreadStruct*  str = (LoadBalancedThreadStruct *)(((itk::MultiThreader::ThreadInfoStruct *)(arg))->UserData);
#endif  

  const int  numberOfEdgesToAnalyze = str->m_Edges.size();
  int  numberOfEdgesAnalyzed = 0;
  while ( true )
    {
    if ( str->m_Builder->GetVerbose() )
      {
      std::cout << "I am thread " << threadId << " and I have analyzed "
                <<  numberOfEdgesAnalyzed << " edges so far" << std::endl;
      }           

    try
      {
      if ( !str->m_Builder->LoadBalancedAnalyzeEdgeFast( str->m_Edges, str->m_PointOccupancies, threadId ) )
        {
        std::cout << "Nothing left to do for thread " << threadId << std::endl;
        break;
        }
      }
    catch( itk::ExceptionObject& e )
      {
      // Apparently somewhere an exception was thrown. We'll just catch it, display what the problem was,
      // and move on
      std::cout << "Exception === Exception === Exception === Exception === Exception === Exception" << std::endl;  
      std::cout << "   An exception was thrown in thread " << threadId << std::endl;
      std::cout << "     " << e << std::endl;
      std::cout << "Exception === Exception === Exception === Exception === Exception === Exception" << std::endl;  
      }

    numberOfEdgesAnalyzed++;
    
    if ( ( threadId == 0 ) &&
         ( numberOfEdgesAnalyzed % 10 == 0 ) )
      {
      //std::cout << "numberOfEdgesAnalyzed: " << numberOfEdgesAnalyzed << std::endl;  
      //std::cout << ( numberOfEdgesAnalyzed % 10 == 0 ) << std::endl;
      
      const double  progress = 1 - static_cast< double >( str->m_Edges.size() ) 
                               / static_cast< double >( numberOfEdgesToAnalyze );
      str->m_Builder->SetProgress( progress );
      str->m_Builder->InvokeEvent( EdgeAnalysisProgressEvent() );
      }  
    }

#if ITK_VERSION_MAJOR >= 5
  return itk::ITK_THREAD_RETURN_DEFAULT_VALUE;
#else    
  return ITK_THREAD_RETURN_VALUE;
#endif  
}



} // end namespace kvl
