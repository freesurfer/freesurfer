#include "kvlAtlasMeshBuilder.h"

#include "kvlAtlasMeshCollectionModelLikelihoodCalculator.h"
#include "kvlAtlasMeshCollectionPositionCostCalculator2.h"
#include "kvlParameterOrderPowellOptimizer.h"
#include "kvlAtlasMeshCollectionFastReferencePositionCost.h"





namespace kvl
{


#if 1
//itk::SimpleFastMutexLock  mutex;
AtlasMeshBuilderMutexLock  mutex;
#endif

//
//
//
AtlasMeshBuilder
::AtlasMeshBuilder()
{
  m_InitialSize[ 0 ] = 10;
  m_InitialSize[ 1 ] = 10;
  m_InitialSize[ 2 ] = 10;
  m_NumberOfUpsamplingSteps = 1;
  m_InitialStiffnesses[ 0 ] = 0.1f;
  m_InitialStiffnesses[ 1 ] = 0.1f;
  m_InitialStiffnesses[ 2 ] = 0.1f;
  m_InitialStiffnesses[ 3 ] = 0.1f;
  m_InitialStiffnesses[ 4 ] = 0.1f;

  m_Mesher = MultiResolutionAtlasMesher::New();

  m_PowellAbsolutePrecision = 1.0f;

  m_NumberOfClasses = 0;

  m_stuckCount = 0;

  // for the collapsed label mapping, we initialize with identity
  for(int i=0; i<256; i++)
  {
    std::vector<unsigned char > v;
    v.push_back((unsigned char) i);
    m_mapCompToComp[i]=v;
  }


  m_Current = 0;

  m_IterationNumber = 0;
  m_MaximumNumberOfIterations = 30;
  m_Progress = 0.0f;
  m_EdgeId = 0;

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
::SetLabelImages( const std::vector< LabelImageType::ConstPointer >& labelImages )
{
  m_LabelImages = labelImages;

  this->SetUp();
}



//
//
//
void
AtlasMeshBuilder
::SetInitialSize( unsigned int* size )
{
  m_InitialSize[ 0 ] = size[ 0 ];
  m_InitialSize[ 1 ] = size[ 1 ];
  m_InitialSize[ 2 ] = size[ 2 ];

  this->SetUp();

}




//
//
//
void
AtlasMeshBuilder
::SetInitialStiffnesses( const float* initialStiffnesses )
{
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
AtlasMeshBuilder
::SetInitialMesh( unsigned int* size, const float* initialStiffnesses )
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
AtlasMeshBuilder
::SetUp()
{

  if ( m_LabelImages.size() == 0 )
  {
    return;
  }

  m_Mesher->SetLabelImages( m_LabelImages );
  m_Mesher->SetNumberOfUpsamplingSteps( m_NumberOfUpsamplingSteps );
  m_Mesher->SetTryToBeSparse( true );
  m_Mesher->SetMapCompToComp(m_mapCompToComp);
  m_Mesher->SetUp( m_InitialSize, m_InitialStiffnesses );

  m_NumberOfClasses = m_Mesher->GetNumberOfClasses();

  std::cout << "AtlasMeshBuilder:  set up with m_NumberOfClasses: " << m_NumberOfClasses << std::endl;

}




//
//
//
const AtlasMeshBuilder::LabelImageType*
AtlasMeshBuilder
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
  AtlasMeshBuilder
  ::Build( AtlasMeshCollection* explicitStartCollection )
  {



    if ( explicitStartCollection == 0 )
    {
      m_Mesher->Go();
      m_Current = const_cast< AtlasMeshCollection* >( m_Mesher->GetCurrentMeshCollection() );
    }
    else
    {
      m_Current = explicitStartCollection;
    }




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


      LoadBalancedThreadStruct  str;
      str.m_Builder = this;
      str.m_Edges = edges;
      str.m_PointOccupancies = pointOccupancies;



      // Do the job
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
  AtlasMeshCollection::Pointer
  AtlasMeshBuilder
  ::ApplyParallellMeshOperations( const AtlasMeshCollection*  collection,
                                  const std::set< AtlasMesh::CellIdentifier >&   disappearingCells,
                                  const std::map< AtlasMesh::PointIdentifier, AtlasMesh::PointIdentifier >&   disappearingPointsLookupTable,
                                  const AtlasMesh::PointsContainer*  newReferencePosition,
                                  const std::vector< AtlasMesh::PointsContainer::Pointer >&  newPositions,
                                  const AtlasMesh::PointDataContainer*  newPointParameters ) const
  {
//    std::cout << "Applying parallell mesh operations...." << std::endl;


    // Copy all positions of the points, except for the second point on the edge, which simply
    // disappears, and first point, which gets the previously determined new position
//    std::cout << "Creating positions" << std::endl;
    std::vector< AtlasMesh::PointsContainer::Pointer >  collapsedPositions;
    AtlasMesh::PointsContainer::Pointer   collapsedReferencePosition = 0;
    for ( unsigned int meshNumber = 0; meshNumber < collection->GetNumberOfMeshes()+1; meshNumber++ )
    {
      AtlasMesh::PointsContainer::ConstPointer  thisPosition;
      AtlasMesh::PointsContainer::ConstPointer  thisNewPosition;
      if ( meshNumber < collection->GetNumberOfMeshes() )
      {
        thisPosition = collection->GetPositions()[ meshNumber ].GetPointer();
        thisNewPosition = newPositions[ meshNumber ].GetPointer();
      }
      else
      {
        thisPosition = collection->GetReferencePosition();
        thisNewPosition = newReferencePosition;
      }

      AtlasMesh::PointsContainer::Pointer  collapsedPosition = AtlasMesh::PointsContainer::New();

      AtlasMesh::PointsContainer::ConstIterator  positionIt = thisPosition->Begin();
      while ( positionIt !=  thisPosition->End() )
      {

        if ( thisNewPosition->IndexExists( positionIt.Index() ) )
        {
          collapsedPosition->InsertElement( positionIt.Index(), thisNewPosition->ElementAt( positionIt.Index() ) );
        }
        else if ( disappearingPointsLookupTable.find( positionIt.Index() ) == disappearingPointsLookupTable.end() )
        {
          collapsedPosition->InsertElement( positionIt.Index(),  positionIt.Value() );
        }

        ++positionIt;
      }


      if ( meshNumber < collection->GetNumberOfMeshes() )
      {
        // Push back
        collapsedPositions.push_back( collapsedPosition );
      }
      else
      {
        collapsedReferencePosition = collapsedPosition;
      }

    }


    // Copy all parameters of the points, except for the second point on the edge, which simply
    // disappears, and the first point on the edge, which will get the previously determined point
    // parameters
//    std::cout << "Creating point parameters" << std::endl;

    AtlasMesh::PointDataContainer::Pointer  collapsedPointParameters = AtlasMesh::PointDataContainer::New();
    AtlasMesh::PointDataContainer::ConstIterator  pointParamIt = collection->GetPointParameters()->Begin();
    while ( pointParamIt != collection->GetPointParameters()->End() )
    {
      if ( newPointParameters->IndexExists( pointParamIt.Index() ) )
      {
        collapsedPointParameters->InsertElement( pointParamIt.Index(),  newPointParameters->ElementAt( pointParamIt.Index() ) );
      }
      else if ( disappearingPointsLookupTable.find( pointParamIt.Index() ) == disappearingPointsLookupTable.end() )
      {
        // Copy
        collapsedPointParameters->InsertElement( pointParamIt.Index(),  pointParamIt.Value() );
      }

      ++pointParamIt;
    }


    // Now loop over all cells, and simply copy unless the cell is disappearing. Also, correct all
    // references to the second edge point (which is disappearing) to a reference to the first edge
    // point
//    std::cout << "Creating cells" << std::endl;
    AtlasMesh::CellsContainer::Pointer  collapsedCells = AtlasMesh::CellsContainer::New();
    for ( AtlasMesh::CellsContainer::ConstIterator cellIt = collection->GetCells()->Begin();
          cellIt != collection->GetCells()->End(); ++cellIt )
    {
      // Check if cell is disappearing
      if ( disappearingCells.find( cellIt.Index() ) != disappearingCells.end() )
      {
        continue;
      }


      const AtlasMesh::CellType*  cell = cellIt.Value();

      // Create a new cell of the correct type
      typedef itk::VertexCell< AtlasMesh::CellType >    VertexCell;
      typedef itk::LineCell< AtlasMesh::CellType >      LineCell;
      typedef itk::TriangleCell< AtlasMesh::CellType >  TriangleCell;
      typedef itk::TetrahedronCell< AtlasMesh::CellType >  TetrahedronCell;

      AtlasMesh::CellAutoPointer newCell;

      if ( cell->GetType() == AtlasMesh::CellType::VERTEX_CELL )
      {
        // Create a new vertex cell
        newCell.TakeOwnership( new VertexCell );
      }
      else if ( cell->GetType() == AtlasMesh::CellType::LINE_CELL  )
      {
        // Create a new line cell
        newCell.TakeOwnership( new LineCell );
      }
      else if ( cell->GetType() == AtlasMesh::CellType::TRIANGLE_CELL  )
      {
        // Create a new triangle cell
        newCell.TakeOwnership( new TriangleCell );
      }
      else
      {
        // Create a new tetrahedron cell
        newCell.TakeOwnership( new TetrahedronCell );
      }


      // Add points to the new cell: Copy in most cases, correct instances of second collapse-edge point
      int localId = 0;
      AtlasMesh::CellType::PointIdConstIterator  pointIt = cell->PointIdsBegin();
      while ( pointIt != cell->PointIdsEnd() )
      {
        AtlasMesh::PointIdentifier  pointId = *pointIt;

        // Correct if necessary
        if ( disappearingPointsLookupTable.find( pointId ) != disappearingPointsLookupTable.end() )
        {
          pointId = disappearingPointsLookupTable.find( pointId )->second;
        }

        //
        newCell->SetPointId( localId, pointId );

        ++localId;
        ++pointIt;
      }


      // Add the cell
      collapsedCells->InsertElement( cellIt.Index(), newCell.ReleaseOwnership() );

    } // End loop over all cells


    // Create a new mesh collection to hold the result.
    AtlasMeshCollection::Pointer  collapsed = AtlasMeshCollection::New();
    collapsed->SetPointParameters( collapsedPointParameters );
    collapsed->SetCells( collapsedCells );
    collapsed->SetReferencePosition( collapsedReferencePosition );
    collapsed->SetPositions( collapsedPositions );
    collapsed->SetK( collection->GetK() );


//    std::cout << "...done!" << std::endl;

    return collapsed;
  }





//
//
//
  float
  AtlasMeshBuilder
  ::GetCurrentDataCost() const
  {
    if ( !m_Current )
    {
      return 0.0f;
    }

    float  dataCost;
    float  dummyAlphasCost;
    this->GetDataCostAndAlphasCost( m_Current, dataCost, dummyAlphasCost );

    return dataCost;
  }


//
//
//
  float
  AtlasMeshBuilder
  ::GetCurrentAlphasCost() const
  {

    if ( !m_Current )
    {
      return 0.0f;
    }

    float  dummyDataCost;
    float  alphasCost;
    this->GetDataCostAndAlphasCost( m_Current, dummyDataCost, alphasCost );

    return alphasCost;
  }



//
//
//
  void
  AtlasMeshBuilder
  ::GetCurrentDataAndAlphasCost( float& currentDataCost, float& currentAlphasCost ) const
  {

    if ( !m_Current )
    {
      currentDataCost = 0.0f;
      currentAlphasCost = 0.0f;

      return;
    }

    this->GetDataCostAndAlphasCost( m_Current, currentDataCost, currentAlphasCost );


  }



//
//
//
  float
  AtlasMeshBuilder
  ::GetCurrentPositionCost() const
  {
    if ( !m_Current )
    {
      return 0.0f;
    }

    return this->GetPositionCost( m_Current );
  }


//
//
//
  float
  AtlasMeshBuilder
  ::GetCurrentCost() const
  {

    if ( !m_Current )
    {
      return 0.0f;
    }

    float  dataCost;
    float  alphasCost;
    this->GetDataCostAndAlphasCost( m_Current, dataCost, alphasCost );
    float  positionCost = this->GetPositionCost( m_Current );

    return dataCost + alphasCost + positionCost;
  }








//
// 
//
  AtlasMeshCollection::Pointer
  AtlasMeshBuilder
  ::TryToCollapseFast( const AtlasMeshCollection* miniCollection,
                       AtlasMesh::CellIdentifier  edgeId,
                       float& miniDataCost, float& miniAlphasCost, float& miniPositionCost,
                       std::set< AtlasMesh::CellIdentifier >&  disappearingCells ) const
  {
    //std::cout << "\n\n\n\n\n==================================================" << std::endl;
    //std::cout << "Trying to collapse edge " << edgeId << std::endl;



    // Get the mesh with the edge collapsed.
    AtlasMesh::CellIdentifier  unifiedVertexId;
    AtlasMeshCollection::Pointer  child;
    if ( !miniCollection->GetCollapsed( edgeId, child, disappearingCells, unifiedVertexId ) )
    {
      return 0;
    }

    const AtlasMesh::PointIdentifier  unifiedPointId =
      *( child->GetCells()->ElementAt( unifiedVertexId )->PointIdsBegin() );


    // Set up the cost calculator
    kvl::AtlasMeshCollectionFastReferencePositionCost::ParametersType  initialPosition( 3 );
    initialPosition[ 0 ] = child->GetReferencePosition()->ElementAt( unifiedPointId )[ 0 ];
    initialPosition[ 1 ] = child->GetReferencePosition()->ElementAt( unifiedPointId )[ 1 ];
    initialPosition[ 2 ] = child->GetReferencePosition()->ElementAt( unifiedPointId )[ 2 ];



#if 1


    // Decide whether or not we're gonna subdivide this hexahedron
    bool  optimizeUnifiedPointReferencePosition = false;

    // Look up the label with the highest alpha in the first point of the miniCollection
    AtlasMesh::PointDataContainer::ConstIterator  pointParamIt = miniCollection->GetPointParameters()->Begin();
    int  maximumAlphaLabelNumber = 0;
    float  maximumAlpha = itk::NumericTraits< float >::min();
    for ( unsigned int classNumber = 0; classNumber < pointParamIt.Value().m_Alphas.Size(); classNumber++ )
    {
      if ( pointParamIt.Value().m_Alphas[ classNumber ] > maximumAlpha )
      {
        maximumAlpha = pointParamIt.Value().m_Alphas[ classNumber ];
        maximumAlphaLabelNumber = classNumber;
      }
    }

    // Look at the alphas in each of the points of miniCollection
    const float threshold = 0.90f;
    for ( ; pointParamIt != miniCollection->GetPointParameters()->End(); ++pointParamIt )
    {
      if ( pointParamIt.Value().m_Alphas[ maximumAlphaLabelNumber ] < threshold )
      {
        optimizeUnifiedPointReferencePosition = true;
      }
    }

#endif

    // Set up cost function
    kvl::AtlasMeshCollectionFastReferencePositionCost::Pointer  costFunction =
      kvl::AtlasMeshCollectionFastReferencePositionCost::New();
    costFunction->SetLabelImages( m_LabelImages );
    costFunction->SetMapCompToComp( (std::vector<unsigned char> *)m_mapCompToComp );
    costFunction->SetPointId( unifiedPointId );
    costFunction->SetInitialMeshCollection( child );
    kvl::AtlasMeshCollectionFastReferencePositionCost::ParametersType  optimalPosition = initialPosition;


#if 1

    if ( !optimizeUnifiedPointReferencePosition )
    {
      // std::cout << "NOT optimizing unified reference position" << std::endl;
    }
    else
    {

#endif


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
        if ( optimizer->GetCurrentCost() != itk::NumericTraits< float >::max() )
        {
          optimalPosition = optimizer->GetCurrentPosition();
          //std::cout << "Changed reference position for the unified point " << unifiedPointId
          //          << " from " << initialPosition << " to " << optimalPosition << std::endl;
        }
        else
        {
          return 0;
        }



      }  // End test if point can move

#if 1
    } // End test if we are going to optimize position
#endif


    // Get the cost at the optimal position
    //std::cout << "Getting cost at optimal position " << optimalPosition << ": ";
    if ( !costFunction->GetValue( optimalPosition, miniDataCost, miniAlphasCost, miniPositionCost ) )
    {
      //std::cout << "not possible" << std::endl;
      return 0;
    }
    //std::cout << "       miniDataCost     : " << miniDataCost << std::endl;
    //std::cout << "       miniAlphasCost   : " << miniAlphasCost << std::endl;
    //std::cout << "       miniPositionCost : " << miniPositionCost << std::endl;


    // Return
    if ( std::isnan( miniDataCost + miniAlphasCost + miniPositionCost ) )
    {
      return 0;
    }
    else
    {
      return const_cast< AtlasMeshCollection* >( costFunction->GetCostCalculationMeshCollection() );
    }

  }









//
//
//
  void
  AtlasMeshBuilder
  ::ExpandCollapse( const AtlasMeshCollection* meshCollection, AtlasMesh::CellIdentifier  edgeId,
                    const AtlasMeshCollection*  optimizedOuterChild,
                    AtlasMeshCollection::Pointer& result, std::set< AtlasMesh::CellIdentifier >*  disappearingCells, AtlasMesh::CellIdentifier& unifiedVertexId ) const
  {
 //   std::cout << "\n\n\n\n\n==================================================" << std::endl;
 //   std::cout << "Trying to expand edge collapse for edge " << edgeId << std::endl;




    // Get the mesh with the edge collapsed.
    if ( !meshCollection->GetCollapsed( edgeId, result, *disappearingCells, unifiedVertexId ) )
    {
      itkExceptionMacro( "Couldn't initialize positions with those estimated in mini-mesh!" );
    }


    // Copy the positions estimated in the mini mesh to the collapsed collection
    std::vector< AtlasMesh::PointsContainer::Pointer >  childPositions =  optimizedOuterChild->GetPositions();
 //   std::cout << "Copying positions from mini mesh collection" << std::endl;
    for ( unsigned int meshNumber = 0; meshNumber < childPositions.size(); meshNumber++ )
    {
      AtlasMesh::PointsContainer::ConstIterator  positionIt = childPositions[ meshNumber ]->Begin();
      while ( positionIt !=  childPositions[ meshNumber ]->End() )
      {
        result->GetPositions()[ meshNumber ]->ElementAt( positionIt.Index() ) =  positionIt.Value();

        ++positionIt;
      }

    }



    // Also copy the reference position from the child to the original mesh after collapsing the edge
    //std::cout << "Copying the reference position from the mini-mesh to the result mesh" << std::endl;
    const AtlasMesh::PointIdentifier  unifiedPointId =
      *( result->GetCells()->ElementAt( unifiedVertexId )->PointIdsBegin() );
    result->GetReferencePosition()->ElementAt( unifiedPointId ) =
      optimizedOuterChild->GetReferencePosition()->ElementAt( unifiedPointId );

    // Force the triangle area parameters to be updated!!!
    result->SetK( result->GetK() );


    // First copy the alphas from the original mesh, and then overwrite the alphas within the child
    for ( AtlasMesh::PointDataContainer::Iterator  paramIt = result->GetPointParameters()->Begin();
          paramIt != result->GetPointParameters()->End();
          ++paramIt )
    {
      paramIt.Value().m_Alphas = meshCollection->GetPointParameters()->ElementAt( paramIt.Index() ).m_Alphas;
    }
    for ( AtlasMesh::PointDataContainer::ConstIterator  paramIt = optimizedOuterChild->GetPointParameters()->Begin();
          paramIt != optimizedOuterChild->GetPointParameters()->End();
          ++paramIt )
    {
      result->GetPointParameters()->ElementAt( paramIt.Index() ).m_Alphas = paramIt.Value().m_Alphas;
    }


  }






//
//
//
  AtlasMeshCollection::Pointer
  AtlasMeshBuilder
  ::TryToRetainFast( const AtlasMeshCollection* miniCollectionConst,
                     AtlasMesh::CellIdentifier  edgeId,
                     float& miniDataCost, float& miniAlphasCost, float& miniPositionCost )
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
      return 0;
    }
    if ( !this->OptimizeReferencePositionFast( miniCollection, vertex1Id, false ) )
    {
      return 0;
    }


    // Measure the cost after optimization
    this->GetDataCostAndAlphasCost( miniCollection, miniDataCost, miniAlphasCost );
    miniPositionCost = this->GetPositionCost( miniCollection );


    // Return
    if ( std::isnan( miniDataCost + miniAlphasCost + miniPositionCost ) )
    {
      return 0;
    }
    else
    {
      return miniCollection;
    }


  }




//
//
//
  void
  AtlasMeshBuilder
  ::ExpandRetain( AtlasMeshCollection* meshCollection, AtlasMesh::CellIdentifier  edgeId,
                  const AtlasMeshCollection*  outerMiniCollection,
                  AtlasMeshCollection::Pointer& result )
  {
//    std::cout << "\n\n\n\n\n==================================================" << std::endl;
//    std::cout << "Trying to expand retained edge " << edgeId << std::endl;

    result = meshCollection;

    // Copy the positions estimated in the mini mesh to the whole collection
    std::vector< AtlasMesh::PointsContainer::Pointer >  miniPositions =  outerMiniCollection->GetPositions();
    //std::cout << "Copying positions from mini mesh collection" << std::endl;
    for ( unsigned int meshNumber = 0; meshNumber < miniPositions.size(); meshNumber++ )
    {
      for ( AtlasMesh::PointsContainer::ConstIterator  positionIt = miniPositions[ meshNumber ]->Begin();
            positionIt !=  miniPositions[ meshNumber ]->End(); ++positionIt )
      {
        result->GetPositions()[ meshNumber ]->ElementAt( positionIt.Index() ) =  positionIt.Value();
      }

    }


    // Retrieve the point ids of the vertices of this edge
    AtlasMesh::CellType::PointIdConstIterator  pointIt =
      outerMiniCollection->GetCells()->ElementAt( edgeId )->PointIdsBegin();
    const AtlasMesh::PointIdentifier  point0Id = *pointIt;
    ++pointIt;
    const AtlasMesh::PointIdentifier  point1Id = *pointIt;


    // Also copy the reference position from the mini mesh to the original mesh
    //std::cout << "Copying the reference positions from the mini-mesh to the result mesh" << std::endl;
    result->GetReferencePosition()->ElementAt( point0Id ) =
      outerMiniCollection->GetReferencePosition()->ElementAt( point0Id );
    result->GetReferencePosition()->ElementAt( point1Id ) =
      outerMiniCollection->GetReferencePosition()->ElementAt( point1Id );

    // Force the triangle area parameters to be updated!!!
    result->SetK( meshCollection->GetK() );


    // Copy the alphas from the mini mesh to the whole mesh
    for ( AtlasMesh::PointDataContainer::ConstIterator  paramIt = outerMiniCollection->GetPointParameters()->Begin();
          paramIt != outerMiniCollection->GetPointParameters()->End();
          ++paramIt )
    {
      result->GetPointParameters()->ElementAt( paramIt.Index() ).m_Alphas = paramIt.Value().m_Alphas;
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
    costFunction->SetLabelImages( m_LabelImages );
    costFunction->SetMapCompToComp( (std::vector<unsigned char> *)m_mapCompToComp );
    costFunction->SetPointId( pointId );
    costFunction->SetInitialMeshCollection( meshCollection );
    kvl::AtlasMeshCollectionFastReferencePositionCost::ParametersType  optimalPosition = initialPosition;

#if 1 /* !!!!!!!!!!!!!!!!! */
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

//      std::cout << "parameterOrder: " << parameterOrder << std::endl;


      // Optimize the reference position
//      std::cout << "Optimizing the position of the unified point " << pointId << std::endl;
      ParameterOrderPowellOptimizer::Pointer  optimizer = ParameterOrderPowellOptimizer::New();
      optimizer->SetCostFunction( costFunction );
      optimizer->SetInitialPosition( initialPosition );
      // optimizer->SetAbsolutePrecisionBrent( m_PowellAbsolutePrecision );
      optimizer->SetStepTolerance( m_PowellAbsolutePrecision );
      // optimizer->SetAbsolutePrecision( m_PowellAbsolutePrecision );
      optimizer->SetParameterOrder( parameterOrder );
      optimizer->StartOptimization();

      // Retrieve the optimal reference position
      if ( optimizer->GetCurrentCost() != itk::NumericTraits< float >::max() )
      {
        optimalPosition = optimizer->GetCurrentPosition();
//        std::cout << "                               Changed reference position for the point " << pointId
//                  << " from " << initialPosition << " to " << optimalPosition << std::endl;
        //const float  initialCost = costFunction->GetValue( initialPosition );
        //const float  optimalCost = costFunction->GetValue( optimalPosition );
        //std::cout << "                                   " << initialCost << "  ->  " << optimalCost << std::endl;
      }
      else
      {
        return false;
      }

    } // End test if we need to optimize
#endif /* !!!!!!!!!!!!!!!!! */


    // Get the optimized outer mini mesh from the costFunction
    if ( costFunction->GetValue( optimalPosition ) == itk::NumericTraits< float >::max() )
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
  unsigned int
  AtlasMeshBuilder
  ::CountNumberOfEdges( const AtlasMeshCollection* meshCollection ) const
  {

    unsigned int numberOfEdges = 0;

    AtlasMesh::CellsContainer::ConstIterator  cellIt = meshCollection->GetCells()->Begin();
    AtlasMesh::CellsContainer::ConstIterator  cellEnd = meshCollection->GetCells()->End();
    while ( cellIt != cellEnd )
    {
      AtlasMesh::CellType*  cell = cellIt.Value();

      if( cell->GetType() == AtlasMesh::CellType::LINE_CELL )
      {
        numberOfEdges++;
      }

      ++cellIt;
    }

    return numberOfEdges;
  }



//
//
//
  float
  AtlasMeshBuilder
  ::GetPositionCost( const AtlasMeshCollection* meshCollection ) const
  {
#if 1
    AtlasMeshCollectionPositionCostCalculator::Pointer  calculator = AtlasMeshCollectionPositionCostCalculator::New();
    calculator->SetMeshCollection( const_cast< AtlasMeshCollection* >( meshCollection ) );
    calculator->SetLabelImages( m_LabelImages );
    calculator->SetMapCompToComp((std::vector<unsigned char> *)m_mapCompToComp); 

#if 0
    if ( this->GetDebug() )
    {
      calculator->DebugOn();
    }
#endif

    return calculator->GetPositionCost();
#else
    return 0.0f;
#endif

  }




//
//
//
  void
  AtlasMeshBuilder
  ::GetDataCostAndAlphasCost( const AtlasMeshCollection* meshCollection, float& dataCost, float& alphasCost ) const
  {

    AtlasMeshCollectionModelLikelihoodCalculator::Pointer  calculator = AtlasMeshCollectionModelLikelihoodCalculator::New();
    calculator->SetMeshCollection( meshCollection );
    calculator->SetLabelImages( m_LabelImages );
    calculator->SetMapCompToComp((std::vector<unsigned char> *)m_mapCompToComp); 

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

#if 1 
    mutex.Lock(); 
/*
    {
      std::ostringstream  descriptionStream;
      descriptionStream << "    [THREAD " << threadId << "] locked mutex because trying to find edge to analyze";
      mutex.DescriptiveLock( descriptionStream.str() );
    }
*/
#endif

    // Check if there is anything left to do
    if ( edges.size() == 0 )
    {
#if 1
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
      AtlasMesh::CellType::PointIdConstIterator  pit
      = m_Current->GetCells()->ElementAt( *it )->PointIdsBegin();
      const AtlasMesh::PointIdentifier  p0Id = *pit;
      ++pit;
      const AtlasMesh::PointIdentifier  p1Id = *pit;

      if ( ( pointOccupancies[ p0Id ] == 0 ) &&
           ( pointOccupancies[ p1Id ] == 0 ) )
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
     std::cout << "    [THREAD " << threadId << "] " << "Encountered non-existing edge " << *it << ". Erasing it" << std::endl;
      edges.erase( *it );
    }


    if ( edgeId == itk::NumericTraits< AtlasMesh::CellIdentifier >::max() )
    {
      // There were edges to be analyzed, but can't work on them at this time cause others or
      // working on it
      std::cout << "    [THREAD " << threadId << "] " << "There are still " << edges.size()
                << " edges but they are all protected at this point" << std::endl;

      m_stuckCount++;
      if (m_stuckCount % 10 == 0)
         std::cout << "    m_stuckCount is at " << m_stuckCount << " (exiting at 20000)" << std::endl; 

      if(m_stuckCount>=20000)
        {
           std::ostringstream  stuckStream;
           stuckStream << "threadsNeverReturning";
           m_Current->Write( stuckStream.str().c_str() );
           exit(-1);
        }

#if 1
      mutex.Unlock();
      //mutex.DescriptiveUnlock();
#endif
      usleep( rand() % 1000000 );
      return true;
    }


    AtlasMesh::CellType::PointIdConstIterator  pit
    = m_Current->GetCells()->ElementAt( edgeId )->PointIdsBegin();
    const AtlasMesh::PointIdentifier  p0Id = *pit;
    ++pit;
    const AtlasMesh::PointIdentifier  p1Id = *pit;
    std::cout << "    [THREAD " << threadId << "] " << "    Analyzing edge with id: " << edgeId
              << " (pointIds " << p0Id << " and " << p1Id << " ) "
              << " (reference positions " << m_Current->GetReferencePosition()->ElementAt( p0Id )
              << " and " << m_Current->GetReferencePosition()->ElementAt( p1Id ) << ")" << std::endl;


    // Remove this edge from the edges to be analyzed
    std::cout << "    [THREAD " << threadId << "] " << "         Removing edge with id: " << edgeId << " from edges (" << edges.size() << ")" << std::endl;

    edges.erase( edgeId );

    if (edges.size() % 1000 == 0)
      std::cout << "    [THREAD " << threadId << "] " << "         Removed edge with id: " << edgeId << " from edges (" << edges.size() << ")" << std::endl;


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

    // Indicate that we are working on them
    for ( std::vector< AtlasMesh::PointIdentifier >::const_iterator  it = affectedPoints.begin();
          it != affectedPoints.end(); ++it )
    {
      //std::cout << "    [THREAD " << threadId << "] " << "        Edge " << edgeId << " is now locking point " << *it << std::endl;
      ( pointOccupancies[ *it ] )++;
    }


#if 1
    mutex.Unlock();
    //mutex.DescriptiveUnlock();
#endif


//    itk::TimeProbe  timeProbe;
//    timeProbe.Start();


    //  Calculate cost when you just optimize this edge's vertex positions
    std::cout << "    [THREAD " << threadId << "] " << "Trying to retain edge " << edgeId << std::endl;
    float  retainedCost = 0;
    float  retainedDataCost = 0;
    float  retainedAlphasCost = 0;
    float  retainedPositionCost = 0;
    AtlasMeshCollection::Pointer  retainedMiniCollection =
      this->TryToRetainFast( miniCollection, edgeId,
                             retainedDataCost, retainedAlphasCost, retainedPositionCost );
    if ( !retainedMiniCollection )
    {
      // Couldn's retain this edge
      retainedCost = itk::NumericTraits< float >::max();
    }
    else
    {
      retainedCost = retainedDataCost + retainedAlphasCost + retainedPositionCost;
    }



    // Calculate the cost of an edge collapse
    std::cout << "    [THREAD " << threadId << "] " << "Trying to collapse edge " << edgeId << std::endl;
    float  collapsedCost = 0;
    float  collapsedDataCost = 0;
    float  collapsedAlphasCost = 0;
    float  collapsedPositionCost = 0;
    std::set< AtlasMesh::CellIdentifier >  collapsedDisappearingCells;
    AtlasMeshCollection::Pointer  collapsedMiniCollection =
      this->TryToCollapseFast( miniCollection, edgeId,
                               collapsedDataCost, collapsedAlphasCost, collapsedPositionCost,
                               collapsedDisappearingCells );
    if ( !collapsedMiniCollection )
    {
      // Couldn't collapse this edge.
      collapsedCost = itk::NumericTraits< float >::max();
    }
    else
    {
      collapsedCost = collapsedDataCost + collapsedAlphasCost + collapsedPositionCost;
    }


    // Evaluate which move is best
    std::cout << "    [THREAD " << threadId << "] " << "                 retainedCost : " << retainedCost
              << "  (" << retainedDataCost << " + " << retainedAlphasCost
              << " + " << retainedPositionCost <<") " << std::endl;
    std::cout << "    [THREAD " << threadId << "] " << "                 collapsedCost: " << collapsedCost
              << "  (" << collapsedDataCost << " + " << collapsedAlphasCost
              << " + " << collapsedPositionCost <<") " << std::endl;
    std::vector< float >  totalCosts;
    totalCosts.push_back( retainedCost );
    totalCosts.push_back( collapsedCost );
    float  minTotalCost = itk::NumericTraits< float >::max();
    int minTotalCostIndex = -1;
    for ( unsigned int i = 0; i < totalCosts.size(); i++ )
    {
      if ( totalCosts[ i ] < minTotalCost )
      {
        minTotalCost = totalCosts[ i ];
        minTotalCostIndex = i;
      }
    }


//    timeProbe.Stop();
//    std::cout << "    [THREAD " << threadId << "] took " << timeProbe.GetMeanTime() << " seconds to evaluate moves" << std::endl;


#if 1
    if ( minTotalCostIndex == -1 )
    {
#if 1
      mutex.Lock();
/*
      {
        std::ostringstream  descriptionStream;
        descriptionStream << "    [THREAD " << threadId << "] locked mutex because impossible configuration encountered";
        mutex.DescriptiveLock( descriptionStream.str() );
      }
*/
#endif

      std::cout << "    [THREAD " << threadId << "] " << "Impossible configuration encountered at eget with id " << edgeId << std::endl;
/*       std::ostringstream  impossibleStream;
      impossibleStream << "impossible_" << edgeId;
      m_Current->Write( impossibleStream.str().c_str() ); */
      //minTotalCostIndex = 0;

      // Now "un-protect" the points we flagged as being worked on
      for ( std::vector< AtlasMesh::PointIdentifier >::const_iterator  it = affectedPoints.begin();
            it != affectedPoints.end(); ++it )
      {
        //std::cout << "    [THREAD " << threadId << "] " << "        Edge " << edgeId << " is now unlocking point " << *it << std::endl;
        ( pointOccupancies[ *it ] )--;
      }

#if 1
      mutex.Unlock();
      //mutex.DescriptiveUnlock();
#endif

      return true;
    }

#endif


#if 1
    mutex.Lock();
/*
    {
      std::ostringstream  descriptionStream;
      descriptionStream << "    [THREAD " << threadId << "] locked mutex because applying best move to the global mesh";
      mutex.DescriptiveLock( descriptionStream.str() );
    }
*/
#endif

    // Do the best move
    if ( minTotalCostIndex == 0 )
    {
      std::cout << "    [THREAD " << threadId << "] " << "        => retaining edge is best solution" << std::endl;

      // this->ExpandRetain( const_cast< AtlasMeshCollection* >( m_Current.GetPointer() ), m_EdgeId,
      //                     m_RetainedMiniCollection,
      //                     newCurrent );


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
      std::cout << "    [THREAD " << threadId << "] " << "        => collapsing edge is best solution" << std::endl;

      // std::set< AtlasMesh::CellIdentifier >  disappearingCells;
      // AtlasMesh::CellIdentifier  unifiedVertexIdDummy;
      // this->ExpandCollapse( m_Current, m_EdgeId, m_CollapsedMiniCollection,
      //                         newCurrent, &disappearingCells, unifiedVertexIdDummy );

      // Look up the ids of the two points on the edge
      AtlasMesh::CellType::PointIdConstIterator  pointIt = m_Current->GetCells()->ElementAt( edgeId )->PointIdsBegin();
      const AtlasMesh::PointIdentifier  edgePoint0Id = *pointIt;
      ++pointIt;
      const AtlasMesh::PointIdentifier  edgePoint1Id = *pointIt;


      // Let's quickly collapse a little bit bigger mesh, and get its cell links. The cell links belonging to the
      // inner points, i.e. the points of our miniCollection, will be copied later on
      AtlasMeshCollection::ConstPointer  biggerMiniCollection = m_Current->GetRegionGrown( edgeId, 2 ).GetPointer();
#if 1
      if ( !biggerMiniCollection )
      {
        std::cout << "    [THREAD " << threadId << "] " << "Ouch!" << std::endl;
        exit( -1 );
      }
#endif
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


    // Now "un-protect" the points we flagged as being worked on
    for ( std::vector< AtlasMesh::PointIdentifier >::const_iterator  it = affectedPoints.begin();
          it != affectedPoints.end(); ++it )
    {
      //std::cout << "    [THREAD " << threadId << "] " << "        Edge " << edgeId << " is now unlocking point " << *it << std::endl;
      ( pointOccupancies[ *it ] )--;
    }

    std::cout << "    [THREAD " << threadId << "] " << "Done with edge " << edgeId << std::endl;

#if 1
    mutex.Unlock();
    //mutex.DescriptiveUnlock();
#endif

    return true;

  }




//
//
//
  void
  AtlasMeshBuilder
  ::AnalyzeEdgeFast( const AtlasMeshCollection*  miniCollection,
                     AtlasMesh::CellIdentifier edgeId,
                     AtlasMesh::PointIdentifier  edgePoint0Id, AtlasMesh::PointIdentifier edgePoint1Id,
                     std::map< AtlasMesh::PointIdentifier, AtlasMesh::PointIdentifier >&  disappearingPointsLookupTable,
                     AtlasMesh::PointsContainer*  newReferencePosition,
                     std::vector< AtlasMesh::PointsContainer::Pointer >& newPositions,
                     AtlasMesh::PointDataContainer*  newPointParameters,
                     std::set< AtlasMesh::CellIdentifier >&  disappearingCells )
  {

//    std::cout << "    Analyzing edge with id: " << edgeId << std::endl;


    //  Calculate cost when you just optimize this edge's vertex positions
    float  retainedCost = 0;
    float  retainedDataCost = 0;
    float  retainedAlphasCost = 0;
    float  retainedPositionCost = 0;
    AtlasMeshCollection::Pointer  retainedMiniCollection =
      this->TryToRetainFast( miniCollection, edgeId,
                             retainedDataCost, retainedAlphasCost, retainedPositionCost );
    if ( !retainedMiniCollection )
    {
      // Couldn's retain this edge
      retainedCost = itk::NumericTraits< float >::max();
    }
    else
    {
      retainedCost = retainedDataCost + retainedAlphasCost + retainedPositionCost;
    }



    // Calculate the cost of an edge collapse
    float  collapsedCost = 0;
    float  collapsedDataCost = 0;
    float  collapsedAlphasCost = 0;
    float  collapsedPositionCost = 0;
    std::set< AtlasMesh::CellIdentifier >  collapsedDisappearingCells;
    AtlasMeshCollection::Pointer  collapsedMiniCollection =
      this->TryToCollapseFast( miniCollection, edgeId,
                               collapsedDataCost, collapsedAlphasCost, collapsedPositionCost,
                               collapsedDisappearingCells );
    if ( !collapsedMiniCollection )
    {
      // Couldn't collapse this edge.
      collapsedCost = itk::NumericTraits< float >::max();
    }
    else
    {
      collapsedCost = collapsedDataCost + collapsedAlphasCost + collapsedPositionCost;
    }


    // Evaluate which move is best
//    std::cout << "                 retainedCost : " << retainedCost
//              << "  (" << retainedDataCost << " + " << retainedAlphasCost
//              << " + " << retainedPositionCost <<") " << std::endl;
 //   std::cout << "                 collapsedCost: " << collapsedCost
//              << "  (" << collapsedDataCost << " + " << collapsedAlphasCost
//              << " + " << collapsedPositionCost <<") " << std::endl;
    std::vector< float >  totalCosts;
    totalCosts.push_back( retainedCost );
    totalCosts.push_back( collapsedCost );
    float  minTotalCost = itk::NumericTraits< float >::max();
    int minTotalCostIndex = -1;
    for ( unsigned int i = 0; i < totalCosts.size(); i++ )
    {
      if ( totalCosts[ i ] < minTotalCost )
      {
        minTotalCost = totalCosts[ i ];
        minTotalCostIndex = i;
      }
    }

#if 1
    if ( minTotalCostIndex == -1 )
    {
#if 1
      mutex.Lock();
#endif

      std::cout << "Impossible configuration encountered at eget with id " << edgeId << std::endl;
      /* std::ostringstream  impossibleStream;
      impossibleStream << "impossible_" << edgeId;
      m_Current->Write( impossibleStream.str().c_str() ); */
      //minTotalCostIndex = 0;

#if 1
      mutex.Unlock();
#endif

      return;
    }

#endif



    // Do the best move
    if ( minTotalCostIndex == 0 )
    {
 //     std::cout << "        => retaining edge is best solution" << std::endl;

      // this->ExpandRetain( const_cast< AtlasMeshCollection* >( m_Current.GetPointer() ), m_EdgeId,
      //                     m_RetainedMiniCollection,
      //                     newCurrent );

      // Remember the reference position of each of the two points
      newReferencePosition->InsertElement( edgePoint0Id,
                                           retainedMiniCollection->GetReferencePosition()->ElementAt( edgePoint0Id ) );
      newReferencePosition->InsertElement( edgePoint1Id,
                                           retainedMiniCollection->GetReferencePosition()->ElementAt( edgePoint1Id ) );

      // Remember the positions of each of the two points
      for ( unsigned int meshNumber = 0; meshNumber < newPositions.size(); meshNumber++ )
      {
        newPositions[ meshNumber ]->InsertElement( edgePoint0Id,
            retainedMiniCollection->GetPositions()[ meshNumber ]->ElementAt( edgePoint0Id ) );
        newPositions[ meshNumber ]->InsertElement( edgePoint1Id,
            retainedMiniCollection->GetPositions()[ meshNumber ]->ElementAt( edgePoint1Id ) );
      }

      // Remember the point parameters of each of the two points
      newPointParameters->InsertElement( edgePoint0Id,
                                         retainedMiniCollection->GetPointParameters()->ElementAt( edgePoint0Id ) );
      newPointParameters->InsertElement( edgePoint1Id,
                                         retainedMiniCollection->GetPointParameters()->ElementAt( edgePoint1Id ) );

    }
    else if ( minTotalCostIndex == 1 )
    {
 //     std::cout << "        => collapsing edge is best solution" << std::endl;

      // std::set< AtlasMesh::CellIdentifier >  disappearingCells;
      // AtlasMesh::CellIdentifier  unifiedVertexIdDummy;
      // this->ExpandCollapse( m_Current, m_EdgeId, m_CollapsedMiniCollection,
      //                         newCurrent, &disappearingCells, unifiedVertexIdDummy );

      // Merge the extra disappearing cells due to this collapse with the cells we already
      // know are disappearing
      std::set< AtlasMesh::CellIdentifier >  totalDisappearingCells;
      std::set_union( disappearingCells.begin(), disappearingCells.end(),
                      collapsedDisappearingCells.begin(), collapsedDisappearingCells.end(),
                      std::inserter( totalDisappearingCells, totalDisappearingCells.begin() ) );
      disappearingCells = totalDisappearingCells;


      // Insert them in the disappearingPointsLookupTable
      disappearingPointsLookupTable[ edgePoint1Id ] = edgePoint0Id;

      // Remember the reference position of the unified point
      newReferencePosition->InsertElement( edgePoint0Id,
                                           collapsedMiniCollection->GetReferencePosition()->ElementAt( edgePoint0Id ) );

      // Remember the positions of the unified point
      for ( unsigned int meshNumber = 0; meshNumber < newPositions.size(); meshNumber++ )
      {
        newPositions[ meshNumber ]->InsertElement( edgePoint0Id,
            collapsedMiniCollection->GetPositions()[ meshNumber ]->ElementAt( edgePoint0Id ) );
      }

      // Remember the point parameters of the unified point
      newPointParameters->InsertElement( edgePoint0Id,
                                         collapsedMiniCollection->GetPointParameters()->ElementAt( edgePoint0Id ) );


    } // End deciding what the best operation is for this edge


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
    std::cout << "Randomizing edge order..." << std::endl;
    std::cout << "   number of edges before: " << edges.size() << std::endl;
    edges =  this->Permutate(  edges );
    std::cout << "   number of edges after: " << edges.size() << std::endl;

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
  std::vector< AtlasMesh::CellIdentifier >
  AtlasMeshBuilder
  ::GetIndependentEdges( int maximumNumberOfIndependentEdges,
                         std::vector< AtlasMesh::CellIdentifier >&  edges ) const
  {

    std::vector< AtlasMesh::CellIdentifier >  edgesToTry;
    if ( edges.size() == 0 )
    {
      edgesToTry = this->GetRandomizedEdges();
    }
    else
    {
      edgesToTry = edges;
    }

    // Select a subset of the edges so that they can all be analyzed independently. Do this
    // by visiting each edge, retrieving a region-grown mesh collection around it, and tagging
    // all points in the region-grown collection as 'untouchable'. When visiting a subsequent
    // edge, this edge will only be added to our subset if it doesn't contain any 'untouchable'
    // points, etc
//    std::cout << "Selecting a subset of edges that can be analyzed independently..." << std::endl;
    std::vector< AtlasMesh::CellIdentifier >  independentEdges;
    std::set< AtlasMesh::PointIdentifier >  untouchablePoints;
    int  numberOfSelectedEdges = 0;
    for ( std::vector< AtlasMesh::CellIdentifier >::const_iterator  it = edgesToTry.begin();
          it != edgesToTry.end(); ++it )
    {
      // Make sure the edge exists still
      if ( !m_Current->GetCells()->IndexExists( *it ) )
      {
        continue;
      }


      // Test first if this edge is independent to previously added edges
      AtlasMesh::CellType::PointIdConstIterator  pit
      = m_Current->GetCells()->ElementAt( *it )->PointIdsBegin();
      const AtlasMesh::PointIdentifier  p0Id = *pit;
      ++pit;
      const AtlasMesh::PointIdentifier  p1Id = *pit;

      if ( ( untouchablePoints.find( p0Id ) != untouchablePoints.end() ) ||
           ( untouchablePoints.find( p1Id ) != untouchablePoints.end() ) )
      {
        // OK, this edge contains an 'untouchable' point. Forget it and move on to
        // the next edge
        continue;
      }


      // Get a region-grown mesh collection around this edge, and tag all its
      // points as 'untouchable'
      AtlasMeshCollection::ConstPointer  regionGrown = m_Current->GetRegionGrown( *it, 1 ).GetPointer();
      for ( AtlasMesh::PointsContainer::ConstIterator  refPosIt = regionGrown->GetReferencePosition()->Begin();
            refPosIt != regionGrown->GetReferencePosition()->End(); ++refPosIt )
      {
        // std::set's insert only inserts if element doesn't exist already, so no worries about double entries
        untouchablePoints.insert( refPosIt.Index() );
      }


      // Remember this edge as an independent one
      independentEdges.push_back( * it );

      //
      numberOfSelectedEdges++;
      if ( numberOfSelectedEdges == maximumNumberOfIndependentEdges )
      {
        break;
      }

    } // End loop over all edges

//    std::cout << "   resulting number of edges: " << independentEdges.size()
//              << " (" << static_cast< float >( independentEdges.size() )
//              / static_cast< float >( edgesToTry.size() ) * 100.0f
 //             << " % of all edges)" << std::endl;


    // Replace the original edges with the set of edges MINUS the ones we returned as independent.
    // We should really replace std::vector with std::set throughout, but I don't wanna break any
    // working code at this point
    std::set< AtlasMesh::CellIdentifier >  edgesToTrySet;
    for ( std::vector< AtlasMesh::CellIdentifier >::const_iterator  it = edgesToTry.begin();
          it != edgesToTry.end(); ++it )
    {
      edgesToTrySet.insert( *it );
    }

    std::set< AtlasMesh::CellIdentifier >  independentEdgesSet;
    for ( std::vector< AtlasMesh::CellIdentifier >::const_iterator  it = independentEdges.begin();
          it != independentEdges.end(); ++it )
    {
      independentEdgesSet.insert( *it );
    }

    std::set< AtlasMesh::CellIdentifier >  untouchedEdgesToTrySet;
    std::set_difference( edgesToTrySet.begin(), edgesToTrySet.end(),
                         independentEdgesSet.begin(), independentEdgesSet.end(),
                         std::inserter( untouchedEdgesToTrySet, untouchedEdgesToTrySet.begin() ) );

    edges.clear();
    for ( std::set< AtlasMesh::CellIdentifier >::const_iterator  it = untouchedEdgesToTrySet.begin();
          it != untouchedEdgesToTrySet.end(); ++it )
    {
      edges.push_back( *it );
    }


    return independentEdges;
  }




//
//
//
  ITK_THREAD_RETURN_TYPE
  AtlasMeshBuilder
  ::LoadBalancedThreaderCallback( void *arg )
  {

    // Retrieve the input arguments
    const int  threadId = ((itk::MultiThreader::ThreadInfoStruct *)(arg))->ThreadID;
    //const int  threadCount = ((itk::MultiThreader::ThreadInfoStruct *)(arg))->NumberOfThreads;

    LoadBalancedThreadStruct*  str = (LoadBalancedThreadStruct *)(((itk::MultiThreader::ThreadInfoStruct *)(arg))->UserData);

    int  numberOfEdgesAnalyzed = 0;
    while ( true )
    {
     std::cout << "I am thread " << threadId << " and I have analyzed "
               <<  numberOfEdgesAnalyzed << " edges so far" << std::endl;

      //if ( !str->m_Builder->LoadBalancedAnalyzeEdge( str->m_Edges, str->m_PointOccupancies, threadId ) )
      if ( !str->m_Builder->LoadBalancedAnalyzeEdgeFast( str->m_Edges, str->m_PointOccupancies, threadId ) )
      {
        std::cout << "Nothing left to do for thread " << threadId << std::endl;
        break;
      }

      numberOfEdgesAnalyzed++;
    }

    return ITK_THREAD_RETURN_VALUE;
  }



} // end namespace kvl
