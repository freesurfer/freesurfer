#include "kvlAtlasMeshBuilder.h"

#include "kvlAtlasMeshCollectionModelLikelihoodCalculator.h"
#include "kvlAtlasMeshCollectionPositionCostCalculator2.h"
#include "kvlParameterOrderPowellOptimizer.h"
#include "kvlAtlasMeshCollectionReferencePositionCost.h"
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



  m_Current = 0;

  m_IterationNumber = 0;
  m_MaximumNumberOfIterations = 30;
  m_Progress = 0.0f;
  m_EdgeId = 0;

  m_RetainedMiniCollection = 0;
  m_CollapsedMiniCollection = 0;
#if 0
  m_SplittedMiniCollection = 0;
  m_SwappedMiniCollection = 0;
#endif

  m_RetainedCost = 0;
  m_RetainedDataCost = 0;
  m_RetainedAlphasCost = 0;
  m_RetainedPositionCost = 0;

  m_CollapsedCost = 0;
  m_CollapsedDataCost = 0;
  m_CollapsedAlphasCost = 0;
  m_CollapsedPositionCost = 0;

#if 0
  m_SplittedCost = 0;
  m_SplittedDataCost = 0;
  m_SplittedAlphasCost = 0;
  m_SplittedPositionCost = 0;

  m_SwappedCost = 0;
  m_SwappedDataCost = 0;
  m_SwappedAlphasCost = 0;
  m_SwappedPositionCost = 0;
#endif


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
#if 0
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
  unsigned int  domainSize[ 3 ];
  domainSize[ 0 ] = static_cast< unsigned int >( labelImageSize[ 0 ] );
  domainSize[ 1 ] = static_cast< unsigned int >( labelImageSize[ 1 ] );
  domainSize[ 2 ] = static_cast< unsigned int >( labelImageSize[ 2 ] );

  m_NumberOfClasses = m_Estimator->GetNumberOfClasses();

  const unsigned int  numberOfMeshes = m_Estimator->GetNumberOfLabelImages();


  // Create a mesh collection accordingly
  AtlasMeshCollection::Pointer  meshCollection = AtlasMeshCollection::New();
  meshCollection->Construct( m_InitialSize, domainSize, m_InitialStiffnesses[ 0 ],
                             m_NumberOfClasses, numberOfMeshes );

  // Now initialize the estimator with the mesh
  m_Estimator->SetInitialMeshCollection( meshCollection );
#else
  if ( m_LabelImages.size() == 0 )
  {
    return;
  }

  m_Mesher->SetLabelImages( m_LabelImages );
  m_Mesher->SetNumberOfUpsamplingSteps( m_NumberOfUpsamplingSteps );
  m_Mesher->SetTryToBeSparse( true );
  m_Mesher->SetUp( m_InitialSize, m_InitialStiffnesses );

  m_NumberOfClasses = m_Mesher->GetNumberOfClasses();

  std::cout << "AtlasMeshBuilder:  set up with m_NumberOfClasses: " << m_NumberOfClasses << std::endl;
#endif

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




#if 0
//
//
//
void
AtlasMeshBuilder
::BuildWithPriorityQueue( AtlasMeshCollection* explicitStartCollection )
{

  if ( explicitStartCollection == 0 )
  {
    // Estimate with original mesh topology, using specified multi-resolution scheme
    AtlasParameterEstimator::LabelImageType::SizeType  labelImageSize =
      m_Estimator->GetLabelImage( 0 )->GetLargestPossibleRegion().GetSize();
    unsigned int  domainSize[ 2 ];
    domainSize[ 0 ] = static_cast< unsigned int >( labelImageSize[ 0 ] );
    domainSize[ 1 ] = static_cast< unsigned int >( labelImageSize[ 1 ] );

    for ( int upsamplingStepNumber = 0; upsamplingStepNumber <= m_NumberOfUpsamplingSteps; upsamplingStepNumber++ )
    {
      // Estimate
      m_Estimator->Estimate();

      // If this is not the final resolution yet, upsample mesh collection
      if ( upsamplingStepNumber != m_NumberOfUpsamplingSteps )
      {
        AtlasMeshCollection::Pointer  upsampledMeshCollection =
          m_Estimator->GetCurrentMeshCollection()->GetUpsampled( domainSize,
              m_InitialStiffnesses[ upsamplingStepNumber + 1 ] );
        m_Estimator->SetInitialMeshCollection( upsampledMeshCollection );

      }

    } // End loop over upsampling steps


  }
  else
  {
    // We have provided an explicit mesh already. Pass it on to the estimator
    m_Estimator->SetInitialMeshCollection( explicitStartCollection );
  }


  // Some initialization stuff
  m_Current = m_Estimator->GetCurrentMeshCollection();
  const unsigned int  maximumPossibleNumberOfEdgeCollapses = this->CountNumberOfEdges( m_Current ) - 5;
  if ( m_MaximumNumberOfEdgeCollapses > maximumPossibleNumberOfEdgeCollapses )
  {
    m_MaximumNumberOfEdgeCollapses = maximumPossibleNumberOfEdgeCollapses;
  }
  m_EdgeCollapseNumber = 0;


#if 0
  for ( unsigned int  reinitializationRoundNumber = 0; ; reinitializationRoundNumber++ )
  {
#if 1
    {
      float beforeDataCost;
      float beforeAlphasCost;
      this->GetDataCostAndAlphasCost( m_Current, beforeDataCost, beforeAlphasCost );
      float beforePositionCost = this->GetPositionCost( m_Current );

      std::cout << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << std::endl;
      std::cout << std::endl;
      std::cout << "Reinitializing at edge collapse number " << m_EdgeCollapseNumber << std::endl;
      std::cout << "    (reinitializationRoundNumber: " << reinitializationRoundNumber << ")" << std::endl;
      std::cout << std::endl;
      std::cout << "Cost just before optimizing for K: " << std::endl;
      std::cout << "       " << beforeDataCost << " + "
      << beforeAlphasCost << " + "
      << beforePositionCost << " = "
      << beforeDataCost + beforeAlphasCost + beforePositionCost << std::endl;
      std::cout << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << std::endl;
    }
#endif

    // Optimize for K
    m_Current = this->OptimizeForK( m_Current );
    //char  dummyChar;
    //std::cout << "Enter a character to continue..." << std::endl;
    //std::cin >> dummyChar;
#if 1
    {
      float afterDataCost;
      float afterAlphasCost;
      this->GetDataCostAndAlphasCost( m_Current, afterDataCost, afterAlphasCost );
      float afterPositionCost = this->GetPositionCost( m_Current );

      std::cout << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << std::endl;
      std::cout << "Cost just after optimizing for K: " << std::endl;
      std::cout << "       " << afterDataCost << " + "
                << afterAlphasCost << " + "
                << afterPositionCost << " = "
                << afterDataCost + afterAlphasCost + afterPositionCost << std::endl;
      std::cout << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << std::endl;
    }
#endif
#endif


    // Initialization
    unsigned int  edgeCollapseNumber = 0;
    AtlasMesh::CellIdentifier   m_LastCollapsedEdgeId = 0;
    AtlasMesh::CellIdentifier  unifiedVertexId = 0;
    m_Gains.clear();
    float  m_CurrentDataCost;
    float  m_CurrentAlphasCost;
    float  m_CurrentPositionCost;
    this->GetDataCostAndAlphasCost( m_Current, m_CurrentDataCost, m_CurrentAlphasCost );
    m_CurrentPositionCost = this->GetPositionCost( m_Current );
    AtlasMeshCollection::ConstPointer  m_Best = m_Current;
    float  m_BestDataCost = m_CurrentDataCost;
    float  m_BestAlphasCost = m_CurrentAlphasCost;
    float  m_BestPositionCost = m_CurrentPositionCost;
    const unsigned int  numberOfEdgeCollapsesToPerform =
      static_cast< unsigned int >( ( this->CountNumberOfEdges( m_Current ) - 5 ) / 1.0f );

#if 0
    if ( reinitializationRoundNumber == 0)
    {
#endif
      this->InvokeEvent( itk::StartEvent() );
#if 0
    }
#endif

    // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    //
    // Main loop: Repeatedly remove one edge
    //
    // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    while ( ( edgeCollapseNumber < numberOfEdgeCollapsesToPerform ) &&
            ( this->CountNumberOfEdges( m_Current ) > 5 ) )
    {


      // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      //
      // Part I. Recalcuate costs of dirty edges
      //
      // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      std::vector< AtlasMesh::CellIdentifier >  edgesToTry;
      if ( ( edgeCollapseNumber == 0 ) || !m_AllowAffectedAreaSpeedUp )
      {
        // The edges that need re-calculation are all the edges
        AtlasMesh::CellsContainer::ConstIterator  cellIt = m_Current->GetCells()->Begin();
        while ( cellIt != m_Current->GetCells()->End() )
        {
          if ( cellIt.Value()->GetType() == AtlasMesh::CellType::LINE_CELL )
          {
            edgesToTry.push_back( cellIt.Index() );
          }

          ++cellIt;
        }

      }
      else
      {
        // The edges that need re-calculation are the edges around the unified vertex
        AtlasMeshCollection::Pointer  affectedCollection =
          m_Current->GetRegionGrown( unifiedVertexId, m_AffectedAreaRadius );
        AtlasMesh::CellsContainer::ConstIterator  cellIt = affectedCollection->GetCells()->Begin();
        while ( cellIt != affectedCollection->GetCells()->End() )
        {
          if ( cellIt.Value()->GetType() == AtlasMesh::CellType::LINE_CELL )
          {
            edgesToTry.push_back( cellIt.Index() );
          }

          ++cellIt;
        }

      }




      // Recalculate gains of all dirty edges
      int  m_NumberOfTriedEdgesInCurrent = 0;
      int  m_NumberOfEdgesToTryInCurrent = edgesToTry.size();
      std::vector< AtlasMesh::CellIdentifier >::const_iterator  edgeIt = edgesToTry.begin();
      while ( edgeIt != edgesToTry.end() )
      {
        // Try the effect of collapsing this edge
        float  dataCost = 0;
        float  alphasCost = 0;
        float  positionCost = 0;
        AtlasMeshCollection::ConstPointer  outerMiniCollection = 0;
        AtlasMeshCollection::ConstPointer  innerMiniCollection = 0;
        if ( !m_AllowMiniMeshSpeedUp )
        {
          // The mini-mesh is the entire mesh
          outerMiniCollection = m_Current;
          innerMiniCollection = m_Current;
        }
        else
        {
          // The mini-mesh is an area surrounding the edge to collapse
          innerMiniCollection = m_Current->GetRegionGrown( *edgeIt, m_MiniMeshRadius );
          outerMiniCollection = m_Current->GetRegionGrown( *edgeIt, m_MiniMeshRadius+1 );
        }

        if ( !this->TryToCollapse( innerMiniCollection, outerMiniCollection,
                                   *edgeIt, dataCost, alphasCost, positionCost ) )
        {
          // remove the edge from the gain container
          m_Gains.Erase( *edgeIt );

          // Prepare for next edge
          m_NumberOfTriedEdgesInCurrent++;
          //this->InvokeEvent( EndEdgeCollapseTrialEvent() );
          ++edgeIt;
          continue;
        }

        float  retainedDataCost;
        float  retainedAlphasCost;
        this->GetDataCostAndAlphasCost( outerMiniCollection, retainedDataCost, retainedAlphasCost );
        const float  retainedPositionCost = this->GetPositionCost( outerMiniCollection );

        // Update the entry of this edge in the gain container
        std::cout << "Setting gain of edge id " << *edgeIt << " to "
                  << retainedDataCost - dataCost << " + " << retainedAlphasCost - alphasCost << " + " << retainedPositionCost - positionCost << " = "
                  <<  retainedDataCost - dataCost + retainedAlphasCost - alphasCost + retainedPositionCost - positionCost << std::endl;
        m_Gains.SetGain( *edgeIt, retainedDataCost - dataCost,
                         retainedAlphasCost - alphasCost,
                         retainedPositionCost - positionCost );


        // Prepare for next edge
        m_NumberOfTriedEdgesInCurrent++;
        //this->InvokeEvent( EndEdgeCollapseTrialEvent() );
        ++edgeIt;
      } // End loop over all dirty edges to update their gains


      //
      //std::cout << "Updated gains container: " << std::endl;
      //m_Gains.Print( std::cout );





      // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      //
      // Part II. Starting from the edge with the highest gain from the gain container, try to collapse it. If an edge won't
      // collapse, remove it. Also remember the id of the unified vertex
      //
      // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      AtlasMesh::CellIdentifier  bestEdgeId = m_Gains.GetBestEdge();
      AtlasMeshCollection::Pointer  newCurrent = AtlasMeshCollection::New();
      std::set< AtlasMesh::CellIdentifier >  disappearingCells;

      while ( true )
      {
        AtlasMeshCollection::ConstPointer  outerMiniCollection = 0;
        AtlasMeshCollection::ConstPointer  innerMiniCollection = 0;
        if ( !m_AllowMiniMeshSpeedUp )
        {
          // The mini-mesh is the entire mesh
          outerMiniCollection = m_Current;
          innerMiniCollection = m_Current;
        }
        else
        {
          // The mini-mesh is an area surrounding the edge to collapse
          innerMiniCollection = m_Current->GetRegionGrown( bestEdgeId, m_MiniMeshRadius );
          outerMiniCollection = m_Current->GetRegionGrown( bestEdgeId, m_MiniMeshRadius+1 );
        }
        float  dataCost = 0;
        float  alphasCost = 0;
        float  positionCost = 0;
        AtlasMeshCollection::Pointer   collapsedMiniCollection =
          this->TryToCollapse( innerMiniCollection, outerMiniCollection,
                               bestEdgeId, dataCost, alphasCost, positionCost );
        if ( !collapsedMiniCollection )
        {
          //
          std::cout << "Best edge in gains container cannot be collapsed. Trying next one (EdgeCollapseNumber: "
                    << m_EdgeCollapseNumber << ")" << std::endl;
          m_Gains.Erase( bestEdgeId );
          if ( m_Gains.size() )
          {
            bestEdgeId = m_Gains.GetBestEdge();
          }
          else
          {
            return;
          }
        }
        else
        {
          this->ExpandCollapse( m_Current, bestEdgeId, collapsedMiniCollection,
                                newCurrent, &disappearingCells, unifiedVertexId );
          break;
        }
      } // End loop over all edges, starting from top of Gains container, until one can be collapsed



      // Replace the current mesh with the new mesh
      std::cout << "Replacing m_Current with newCurrent " << newCurrent.GetPointer() << std::endl;
      m_Current = newCurrent;
      const float  previousDataCost = m_CurrentDataCost;
      const float  previousAlphasCost = m_CurrentAlphasCost;
      const float  previousPositionCost = m_CurrentPositionCost;
      this->GetDataCostAndAlphasCost( m_Current, m_CurrentDataCost, m_CurrentAlphasCost );
      m_CurrentPositionCost = this->GetPositionCost( m_Current );

      std::ostringstream  currentFileNameStream;
      currentFileNameStream << "Current_" << edgeCollapseNumber << "_"  << m_CurrentDataCost << "_" << m_CurrentAlphasCost << "_"
                            << m_CurrentPositionCost << "_____" << m_CurrentDataCost + m_CurrentAlphasCost + m_CurrentPositionCost << ".txt";
      std::cout << "Writing out to file " << currentFileNameStream.str() << std::endl;
      m_Current->Write( currentFileNameStream.str().c_str() );

#if 0
      //
      std::ostringstream  debugConsistencyStream;
      debugConsistencyStream << "debugConsistency_iteration" <<  m_EdgeCollapseNumber
                             << "_edge" << bestEdgeId;

      std::ofstream  out( debugConsistencyStream.str().c_str(), std::ios::app );
      if ( !out.bad() )
      {
        out << "Collapsing edge with id: " << bestEdgeId << std::endl;
        float  promisedDataGain;
        float  promisedAlphasGain;
        float  promisedPositionGain;
        m_Gains.GetGain( bestEdgeId, promisedDataGain, promisedAlphasGain, promisedPositionGain );
        out << "      Gain container predicted: "
            << promisedDataGain << " " << promisedAlphasGain << " " << promisedPositionGain
            << " = " << promisedDataGain + promisedAlphasGain + promisedPositionGain << std::endl;
        out << "      We actually got: "
            << dataGain << " " << alphasGain << " " << positionGain
            << " = " << dataGain + alphasGain + positionGain  << std::endl;
        out << "      The final end result is : "
            << previousDataCost - m_CurrentDataCost << " "
            << previousAlphasCost - m_CurrentAlphasCost << " "
            << previousPositionCost - m_CurrentPositionCost
            << " = " <<  previousDataCost + previousAlphasCost + previousPositionCost
            - m_CurrentDataCost - m_CurrentAlphasCost - m_CurrentPositionCost << std::endl;
      }
#endif


#if 0
      // Do a consistency check, comparing the minimesh data gain the gain container promised, the minimesh
      // gain we actually got, and the gain calculated on the whole mesh
      float  promisedDataGain;
      float  promisedAlphasGain;
      float  promisedPositionGain;
      m_Gains.GetGain( bestEdgeId, promisedDataGain, promisedAlphasGain, promisedPositionGain );
      const float  promisedGain = promisedDataGain + promisedAlphasGain + promisedPositionGain;
      const float  actualGain = dataGain + alphasGain + positionGain;
      const float  finalEndGain = previousDataCost + previousAlphasCost + previousPositionCost
                                  - m_CurrentDataCost - m_CurrentAlphasCost - m_CurrentPositionCost;
      const float  averageGain = ( promisedGain + actualGain + finalEndGain ) / 3;
      const float  promisedGainDeviation = ( averageGain - promisedGain ) / averageGain;
      const float  actualGainDeviation = ( averageGain - actualGain ) / averageGain;
      const float  finalEndGainDeviation = ( averageGain - finalEndGain ) / averageGain;
      const float  deviationTolerance = .05;
      if ( ( vnl_math_abs( promisedGainDeviation ) > deviationTolerance ) ||
           ( vnl_math_abs( actualGainDeviation ) > deviationTolerance ) ||
           ( vnl_math_abs( finalEndGainDeviation ) > deviationTolerance ) )
      {
        std::cerr << "Found inconsistency at iteration "
                  << m_EdgeCollapseNumber
                  << " with collapsing edge " << bestEdgeId  << "\n"
                  << "      Gain container predicted: "
                  << promisedDataGain << " " << promisedAlphasGain << " " << promisedPositionGain
                  << " = " << promisedGain << "\n"
                  << "      We actually got: "
                  << dataGain << " " << alphasGain << " " << positionGain
                  << " = " << actualGain  << "\n"
                  << "      The final end result is : "
                  << previousDataCost - m_CurrentDataCost << " "
                  << previousAlphasCost - m_CurrentAlphasCost << " "
                  << previousPositionCost - m_CurrentPositionCost
                  << " = " <<  finalEndGain << std::endl;
      }

#endif


      // remove the cells that have disappeared from the gain container
      std::cout << "Removing cells that have disappeard from the gain container" << std::endl;
      std::set< AtlasMesh::CellIdentifier >::const_iterator  disappearIt = disappearingCells.begin();
      while ( disappearIt != disappearingCells.end() )
      {
        m_Gains.Erase( *disappearIt );
        ++disappearIt;
      }

      // Remember which edge we just collapsed
      m_LastCollapsedEdgeId = bestEdgeId;

      // If the cost is the overall best so far, remember it.
      if ( ( m_CurrentDataCost + m_CurrentAlphasCost + m_CurrentPositionCost ) <
           ( m_BestDataCost + m_BestAlphasCost + m_BestPositionCost ) )
      {
        //m_Best = m_Current;
        m_BestDataCost = m_CurrentDataCost;
        m_BestAlphasCost = m_CurrentAlphasCost;
        m_BestPositionCost = m_CurrentPositionCost;
      }

      // Next iteration
      edgeCollapseNumber++;
      m_EdgeCollapseNumber++;
      //this->InvokeEvent( EdgeCollapseEvent() );
    }

    this->InvokeEvent( itk::EndEvent() );


  }
#endif





#if 0
//
//
//
  void
  AtlasMeshBuilder
  ::Build( AtlasMeshCollection* explicitStartCollection )
  {

    if ( explicitStartCollection == 0 )
    {
      // Estimate with original mesh topology, using specified multi-resolution scheme
      AtlasParameterEstimator::LabelImageType::SizeType  labelImageSize =
        m_Estimator->GetLabelImage( 0 )->GetLargestPossibleRegion().GetSize();
      unsigned int  domainSize[ 3 ];
      domainSize[ 0 ] = static_cast< unsigned int >( labelImageSize[ 0 ] );
      domainSize[ 1 ] = static_cast< unsigned int >( labelImageSize[ 1 ] );
      domainSize[ 2 ] = static_cast< unsigned int >( labelImageSize[ 2 ] );

      for ( int upsamplingStepNumber = 0; upsamplingStepNumber <= m_NumberOfUpsamplingSteps; upsamplingStepNumber++ )
      {
        // Estimate
        m_Estimator->Estimate();

        // If this is not the final resolution yet, upsample mesh collection
        if ( upsamplingStepNumber != m_NumberOfUpsamplingSteps )
        {
          AtlasMeshCollection::Pointer  upsampledMeshCollection =
            m_Estimator->GetCurrentMeshCollection()->GetUpsampled();
          upsampledMeshCollection->SetK( m_InitialStiffnesses[ upsamplingStepNumber + 1 ] );
          m_Estimator->SetInitialMeshCollection( upsampledMeshCollection );

        }

      } // End loop over upsampling steps


    }
    else
    {
      // We have provided an explicit mesh already. Pass it on to the estimator
      m_Estimator->SetInitialMeshCollection( explicitStartCollection );
    }


    // Some initialization stuff
    m_Current = m_Estimator->GetCurrentMeshCollection();

    // ICM-like algorithm
    this->InvokeEvent( itk::StartEvent() );
    for ( m_IterationNumber = 0; m_IterationNumber < m_MaximumNumberOfIterations; m_IterationNumber++ )
    {
      this->InvokeEvent( itk::IterationEvent() );

      std::cout << "Iteration " << m_IterationNumber << std::endl;
#if 0
      {
        std::ostringstream  outFileNameStream;
        outFileNameStream << "builderMeshCollectionAtIteration" << m_IterationNumber << ".txt";
        m_Current->Write( outFileNameStream.str().c_str() );
      }
#endif

      //m_Current = this->OptimizeForK( m_Current );


      //
      // Part I: optimize connectivity
      //
      // Cache edges at this point
      std::vector< AtlasMesh::CellIdentifier >  edgesToTry;
      for ( AtlasMesh::CellsContainer::ConstIterator  cellIt = m_Current->GetCells()->Begin();
            cellIt != m_Current->GetCells()->End(); ++cellIt )
      {
        if ( cellIt.Value()->GetType() == AtlasMesh::CellType::LINE_CELL )
        {
          edgesToTry.push_back( cellIt.Index() );
        }
      }

      // Randomize the edges' order
      std::cout << "Randomizing edge order..." << std::endl;
      std::cout << "   number of edges before: " << edgesToTry.size() << std::endl;
      edgesToTry =  this->Permutate(  edgesToTry );
      std::cout << "   number of edges after: " << edgesToTry.size() << std::endl;

      // Loop over all edges
      unsigned int  progressCounter = 0;
      m_Progress = 0.0f;
      for ( std::vector< AtlasMesh::CellIdentifier >::const_iterator  edgeIt = edgesToTry.begin();
            edgeIt != edgesToTry.end(); ++edgeIt, ++progressCounter )
      {
        m_EdgeId = *edgeIt;
        m_Progress = static_cast< float >( progressCounter ) / static_cast< float >( edgesToTry.size() );


        // Check if the edge is still present in the mesh.
        if ( !m_Current->GetCells()->IndexExists( m_EdgeId ) )
        {
          continue;
        }


        std::cout << "    Analyzing edge with id: " << m_EdgeId << std::endl;
        this->InvokeEvent( StartEdgeAnalysisEvent() );

#if 0
        {
          m_Current->Write( "m_Current.txt" );
        }
#endif


        // Get the mini collection surrounding the edge to try to collapse
        AtlasMeshCollection::ConstPointer  outerMiniCollection = m_Current->GetRegionGrown( m_EdgeId, 2 ).GetPointer();
        AtlasMeshCollection::ConstPointer  innerMiniCollection = m_Current->GetRegionGrown( m_EdgeId, 1 ).GetPointer();


        //  Calculate cost when you just optimize this edge's vertex positions
        m_RetainedCost = 0;
        m_RetainedDataCost = 0;
        m_RetainedAlphasCost = 0;
        m_RetainedPositionCost = 0;
        m_RetainedMiniCollection =
          this->TryToRetain( innerMiniCollection, outerMiniCollection, m_EdgeId,
                             m_RetainedDataCost, m_RetainedAlphasCost, m_RetainedPositionCost );
        if ( !m_RetainedMiniCollection )
        {
          // Couldn's retain this edge
          m_RetainedCost = itk::NumericTraits< float >::max();
        }
        else
        {
          m_RetainedCost = m_RetainedDataCost + m_RetainedAlphasCost + m_RetainedPositionCost;
        }



        // Calculate the cost of an edge collapse
        m_CollapsedCost = 0;
        m_CollapsedDataCost = 0;
        m_CollapsedAlphasCost = 0;
        m_CollapsedPositionCost = 0;
#if 1 /* !!!!!!!!!!!!!!!!!!!!!! */
        m_CollapsedMiniCollection =
          this->TryToCollapse( innerMiniCollection, outerMiniCollection, m_EdgeId,
                               m_CollapsedDataCost, m_CollapsedAlphasCost, m_CollapsedPositionCost );
        if ( !m_CollapsedMiniCollection )
        {
          // Couldn't collapse this edge.
          m_CollapsedCost = itk::NumericTraits< float >::max();
        }
        else
        {
          m_CollapsedCost = m_CollapsedDataCost + m_CollapsedAlphasCost + m_CollapsedPositionCost;
        }
#else /* !!!!!!!!!!!!!!!!!!!!!! */
        m_CollapsedCost = itk::NumericTraits< float >::max();
#endif  /* !!!!!!!!!!!!!!!!!!!!!! */

#if 0
        // Calculate cost of an edge split
        m_SplittedCost = 0;
        m_SplittedDataCost = 0;
        m_SplittedAlphasCost = 0;
        m_SplittedPositionCost = 0;
#if 0 /* !!!!!!!!!!!!!!!!!!!!!! */
        m_SplittedMiniCollection =
          this->TryToSplit( m_Current, m_EdgeId,
                            m_SplittedDataCost, m_SplittedAlphasCost, m_SplittedPositionCost,
                            dummyResult );
        if ( !m_SplittedMiniCollection )
        {
          // Couldn's split this edge
          m_SplittedCost = itk::NumericTraits< float >::max();
        }
        else
        {
          m_SplittedCost = m_SplittedDataCost + m_SplittedAlphasCost + m_SplittedPositionCost;
        }
#else /* !!!!!!!!!!!!!!!!!!!!!! */
        m_SplittedCost = itk::NumericTraits< float >::max();
#endif /* !!!!!!!!!!!!!!!!!!!!!! */


        // Calculate cost of an edge swap
        m_SwappedCost = 0;
        m_SwappedDataCost = 0;
        m_SwappedAlphasCost = 0;
        m_SwappedPositionCost = 0;
#if 0 /* !!!!!!!!!!!!!!!!!!!!!! */
        m_SwappedMiniCollection =
          this->TryToSwap( m_Current, m_EdgeId,
                           m_SwappedDataCost, m_SwappedAlphasCost, m_SwappedPositionCost,
                           dummyResult );
        if ( !m_SwappedMiniCollection )
        {
          // Couldn's swap this edge
          m_SwappedCost = itk::NumericTraits< float >::max();
        }
        else
        {
          m_SwappedCost = m_SwappedDataCost + m_SwappedAlphasCost + m_SwappedPositionCost;
        }
#else /* !!!!!!!!!!!!!!!!!!!!!! */
        m_SwappedCost = itk::NumericTraits< float >::max();
#endif /* !!!!!!!!!!!!!!!!!!!!!! */
#endif

        // Evaluate which move is best
        std::cout << "                 retainedCost : " << m_RetainedCost
                  << "  (" << m_RetainedDataCost << " + " << m_RetainedAlphasCost
                  << " + " << m_RetainedPositionCost <<") " << std::endl;
        std::cout << "                 collapsedCost: " << m_CollapsedCost
                  << "  (" << m_CollapsedDataCost << " + " << m_CollapsedAlphasCost
                  << " + " << m_CollapsedPositionCost <<") " << std::endl;
#if 0
        std::cout << "                 splittedCost : " << m_SplittedCost
                  << "  (" << m_SplittedDataCost << " + " << m_SplittedAlphasCost
                  << " + " << m_SplittedPositionCost <<") " << std::endl;
        std::cout << "                 swappedCost  : " << m_SwappedCost
                  << "  (" << m_SwappedDataCost << " + " << m_SwappedAlphasCost
                  << " + " << m_SwappedPositionCost <<") " << std::endl;
#endif
        std::vector< float >  totalCosts;
        totalCosts.push_back( m_RetainedCost );
        totalCosts.push_back( m_CollapsedCost );
#if 0
        totalCosts.push_back( m_SplittedCost );
        totalCosts.push_back( m_SwappedCost );
#endif
        float  minTotalCost = itk::NumericTraits< float >::max();
        int minTotalCostIndex = -1;
        for ( int i = 0; i < totalCosts.size(); i++ )
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
          std::cout << "Impossible configuration encountered at eget with id " << m_EdgeId << std::endl;
          std::ostringstream  impossibleStream;
          impossibleStream << "impossible_" << m_EdgeId;
          m_Current->Write( impossibleStream.str().c_str() );
          //minTotalCostIndex = 0;
          continue;
        }

#endif


        // Do the best move
        AtlasMeshCollection::Pointer  newCurrent = AtlasMeshCollection::New();
        if ( minTotalCostIndex == 0 )
        {
          std::cout << "        => retaining edge is best solution" << std::endl;
          this->ExpandRetain( const_cast< AtlasMeshCollection* >( m_Current.GetPointer() ), m_EdgeId,
                              m_RetainedMiniCollection,
                              newCurrent );
        }
        else if ( minTotalCostIndex == 1 )
        {
          std::cout << "        => collapsing edge is best solution" << std::endl;

          std::set< AtlasMesh::CellIdentifier >  disappearingCells;
          AtlasMesh::CellIdentifier  unifiedVertexIdDummy;
          this->ExpandCollapse( m_Current, m_EdgeId, m_CollapsedMiniCollection,
                                newCurrent, &disappearingCells, unifiedVertexIdDummy );
        }
#if 0
        else if ( minTotalCostIndex == 2 )
        {
          std::cout << "        => splitting edge is best solution" << std::endl;

          this->TryToSplit( m_Current, m_EdgeId,
                            m_SplittedDataCost, m_SplittedAlphasCost, m_SplittedPositionCost,
                            newCurrent );
        }
        else
        {
          std::cout << "        => swapping edge is best solution" << std::endl;

          this->TryToSwap( m_Current, m_EdgeId,
                           m_SwappedDataCost, m_SwappedAlphasCost, m_SwappedPositionCost,
                           newCurrent );
        }
#endif







        // Replace the current mesh with the new mesh
        //std::cout << "Replacing m_Current with newCurrent " << newCurrent.GetPointer() << std::endl;
        m_Current = newCurrent;

#if 1
        {
          std::ostringstream  outFileNameStream;
          outFileNameStream << "builderMeshCollection_Iteration" << m_IterationNumber << "_edge" << m_EdgeId << ".txt";
          m_Current->Write( outFileNameStream.str().c_str() );
        }
#endif

        // Next edge
        this->InvokeEvent( EndEdgeAnalysisEvent() );

      } // End loop over all edges



#if 0
      //
      // Re-estimate mesh again from start position
      //
      std::cout << "Re-estimating positions and alpahs from scratch" << std::endl;

      // Assign flat priors
      AtlasMeshCollection::Pointer  nonConstCurrent = const_cast< AtlasMeshCollection* >( m_Current.GetPointer() );
      AtlasMesh::PointDataContainer::Iterator  pointParamIt = nonConstCurrent->GetPointParameters()->Begin();
      const unsigned int  numberOfClasses = pointParamIt.Value().m_Alphas.Size();
      AtlasAlphasType   flatAlphasEntry( numberOfClasses );
      flatAlphasEntry.Fill( 1.0f / static_cast< float >( numberOfClasses ) );
      for ( ; pointParamIt != nonConstCurrent->GetPointParameters()->End(); ++pointParamIt )
      {
        if ( pointParamIt.Value().m_CanChangeAlphas )
        {
          pointParamIt.Value().m_Alphas = flatAlphasEntry;
        }
      }

      // Set positions to reference position
      for ( unsigned int  meshNumber=0; meshNumber < nonConstCurrent->GetNumberOfMeshes(); meshNumber++ )
      {
        AtlasMesh::PointsContainer::ConstIterator  refPosIt = nonConstCurrent->GetReferencePosition()->Begin();
        AtlasMesh::PointsContainer::Iterator  posIt = nonConstCurrent->GetPositions()[ meshNumber ]->Begin();
        for ( ; refPosIt != nonConstCurrent->GetReferencePosition()->End(); ++refPosIt, ++posIt )
        {
          posIt.Value() = refPosIt.Value();
        }

      }

      // Rebuild cache
      nonConstCurrent->SetK( nonConstCurrent->GetK() );

      // Estimate
      m_Estimator->SetInitialMeshCollection( nonConstCurrent );
      m_Estimator->Estimate();
      m_Current = m_Estimator->GetCurrentMeshCollection();
#endif


      //
      // Part II: optimize geometry
      //
#if 0
      // Loop over all vertices
      for ( AtlasMesh::CellsContainer::ConstIterator  cellIt = m_Current->GetCells()->Begin();
            cellIt != m_Current->GetCells()->End(); ++cellIt )
      {
        if ( cellIt.Value()->GetType() != AtlasMesh::CellType::VERTEX_CELL )
        {
          continue;
        }

        // Optimize
        std::cout << "    Optimizing reference position of vertex with id: " << cellIt.Index() << std::endl;
        /*      {
              float dataCost;
              float alphasCost;
              this->GetDataCostAndAlphasCost( m_Current, dataCost, alphasCost );
              float positionCost = this->GetPositionCost( m_Current );
              std::cout << "                 cost before: " << dataCost + alphasCost + positionCost << std::endl;
              m_Current->Write( "before.txt" );
              }*/
        this->OptimizeReferencePosition( const_cast< AtlasMeshCollection* >( m_Current.GetPointer() ), cellIt.Index() );
//       {
//       float dataCost;
//       float alphasCost;
//       this->GetDataCostAndAlphasCost( m_Current, dataCost, alphasCost );
//       float positionCost = this->GetPositionCost( m_Current );
//       std::cout << "                 cost after: " << dataCost + alphasCost + positionCost << std::endl;
//       m_Current->Write( "after.txt" );
//       }

        this->InvokeEvent( EdgeCollapseEvent() );
      } // End loop over vertices
#endif



#if 0 /* !!!!!!!!!!!!!!!!!!!!!! */
      //
      // Part II: optimize K
      //
      m_Current = this->OptimizeForK( m_Current );
#endif /* !!!!!!!!!!!!!!!!!!!!!! */

    } // End ICM iterations

    this->InvokeEvent( itk::EndEvent() );


  }

#elif 0


//
//
//
  void
  AtlasMeshBuilder
  ::Build( AtlasMeshCollection* explicitStartCollection )
  {

    if ( explicitStartCollection == 0 )
    {
      // Estimate with original mesh topology, using specified multi-resolution scheme
      AtlasParameterEstimator::LabelImageType::SizeType  labelImageSize =
        m_Estimator->GetLabelImage( 0 )->GetLargestPossibleRegion().GetSize();
      unsigned int  domainSize[ 3 ];
      domainSize[ 0 ] = static_cast< unsigned int >( labelImageSize[ 0 ] );
      domainSize[ 1 ] = static_cast< unsigned int >( labelImageSize[ 1 ] );
      domainSize[ 2 ] = static_cast< unsigned int >( labelImageSize[ 2 ] );

      for ( int upsamplingStepNumber = 0; upsamplingStepNumber <= m_NumberOfUpsamplingSteps; upsamplingStepNumber++ )
      {
        // Estimate
        m_Estimator->Estimate();

        // If this is not the final resolution yet, upsample mesh collection
        if ( upsamplingStepNumber != m_NumberOfUpsamplingSteps )
        {
          AtlasMeshCollection::Pointer  upsampledMeshCollection =
            m_Estimator->GetCurrentMeshCollection()->GetUpsampled();
          upsampledMeshCollection->SetK( m_InitialStiffnesses[ upsamplingStepNumber + 1 ] );
          m_Estimator->SetInitialMeshCollection( upsampledMeshCollection );

        }

      } // End loop over upsampling steps


    }
    else
    {
      // We have provided an explicit mesh already. Pass it on to the estimator
      m_Estimator->SetInitialMeshCollection( explicitStartCollection );
    }


    // Some initialization stuff
    m_Current = m_Estimator->GetCurrentMeshCollection();

#if 1
    std::vector< AtlasMesh::CellIdentifier >  debugPreviousEdgesToTry;
#endif

    // ICM-like algorithm
    this->InvokeEvent( itk::StartEvent() );
    for ( m_IterationNumber = 0; m_IterationNumber < m_MaximumNumberOfIterations; m_IterationNumber++ )
    {
      this->InvokeEvent( itk::IterationEvent() );

      std::cout << "Iteration " << m_IterationNumber << std::endl;

      // Cache edges at this point
      std::vector< AtlasMesh::CellIdentifier >  edgesToTry;
      for ( AtlasMesh::CellsContainer::ConstIterator  cellIt = m_Current->GetCells()->Begin();
            cellIt != m_Current->GetCells()->End(); ++cellIt )
      {
        if ( cellIt.Value()->GetType() == AtlasMesh::CellType::LINE_CELL )
        {
          edgesToTry.push_back( cellIt.Index() );
        }
      }

      // Randomize the edges' order
      std::cout << "Randomizing edge order..." << std::endl;
      std::cout << "   number of edges before: " << edgesToTry.size() << std::endl;
      edgesToTry =  this->Permutate(  edgesToTry );
      std::cout << "   number of edges after: " << edgesToTry.size() << std::endl;


      // Select a subset of the edges so that they can all be analyzed independently. Do this
      // by visiting each edge, retrieving a region-grown mesh collection around it, and tagging
      // all points in the region-grown collection as 'untouchable'. When visiting a subsequent
      // edge, this edge will only be added to our subset if it doesn't contain any 'untouchable'
      // points, etc
      std::cout << "Selecting a subset of edges that can be analyzed independently..." << std::endl;
      std::vector< AtlasMesh::CellIdentifier >  independentEdges;
      std::set< AtlasMesh::PointIdentifier >  untouchablePoints;
#if 1
      int  numberOfSelectedEdges = 0;
#endif
      for ( std::vector< AtlasMesh::CellIdentifier >::const_iterator  it = edgesToTry.begin();
            it != edgesToTry.end(); ++it )
      {
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

#if 1
        //
        numberOfSelectedEdges++;
        if ( numberOfSelectedEdges == 20 )
        {
          break;
        }
#endif

      } // End loop over all edges

      std::cout << "   resulting number of edges: " << independentEdges.size()
                << " (" << static_cast< float >( independentEdges.size() )
                / static_cast< float >( edgesToTry.size() ) * 100.0f
                << " % of all edges)" << std::endl;
      edgesToTry = independentEdges;


#if 1
      int  numberOfReoccuringEdges = 0;
      std::set< AtlasMesh::CellIdentifier >  debugSet;
      for ( std::vector< AtlasMesh::CellIdentifier >::const_iterator  edgeIt = debugPreviousEdgesToTry.begin();
            edgeIt != debugPreviousEdgesToTry.end(); ++edgeIt )
      {
        debugSet.insert( *edgeIt );
      }

      for ( std::vector< AtlasMesh::CellIdentifier >::const_iterator  edgeIt = edgesToTry.begin();
            edgeIt != edgesToTry.end(); ++edgeIt )
      {
        if ( debugSet.find( *edgeIt ) != debugSet.end() )
        {
          numberOfReoccuringEdges++;
        }
      }
      std::cout <<  static_cast< float >( numberOfReoccuringEdges )
                / static_cast< float >( edgesToTry.size() ) * 100.0f
                << "% of these edges were also visited last iteration"<< std::endl;

      debugPreviousEdgesToTry = edgesToTry;
#endif


      // Loop over all edges
      unsigned int  progressCounter = 0;
      m_Progress = 0.0f;
      std::map< AtlasMesh::PointIdentifier, AtlasMesh::PointIdentifier >   disappearingPointsLookupTable;
      AtlasMesh::PointsContainer::Pointer   newReferencePosition = AtlasMesh::PointsContainer::New();
      std::vector< AtlasMesh::PointsContainer::Pointer >   newPositions;
      int  numberOfEdgeCollapses = 0;
      for ( unsigned int meshNumber = 0; meshNumber < m_Current->GetNumberOfMeshes(); meshNumber++ )
      {
        newPositions.push_back( AtlasMesh::PointsContainer::New() );
      }
      AtlasMesh::PointDataContainer::Pointer  newPointParameters = AtlasMesh::PointDataContainer::New();
      std::set< AtlasMesh::CellIdentifier >   disappearingCells;
      for ( std::vector< AtlasMesh::CellIdentifier >::const_iterator  edgeIt = edgesToTry.begin();
            edgeIt != edgesToTry.end(); ++edgeIt, ++progressCounter )
      {
        m_EdgeId = *edgeIt;
        m_Progress = static_cast< float >( progressCounter ) / static_cast< float >( edgesToTry.size() );


        // // Check if the edge is still present in the mesh. --> This shouldn't be necessary!!!
        // if ( !m_Current->GetCells()->IndexExists( m_EdgeId ) )
        //   {
        //   continue;
        //   }


        std::cout << "    Analyzing edge with id: " << m_EdgeId << std::endl;
        this->InvokeEvent( StartEdgeAnalysisEvent() );


        // Get the mini collection surrounding the edge to try to collapse
        AtlasMeshCollection::ConstPointer  outerMiniCollection = m_Current->GetRegionGrown( m_EdgeId, 2 ).GetPointer();
        AtlasMeshCollection::ConstPointer  innerMiniCollection = m_Current->GetRegionGrown( m_EdgeId, 1 ).GetPointer();


        //  Calculate cost when you just optimize this edge's vertex positions
        m_RetainedCost = 0;
        m_RetainedDataCost = 0;
        m_RetainedAlphasCost = 0;
        m_RetainedPositionCost = 0;
        m_RetainedMiniCollection =
          this->TryToRetain( innerMiniCollection, outerMiniCollection, m_EdgeId,
                             m_RetainedDataCost, m_RetainedAlphasCost, m_RetainedPositionCost );
        if ( !m_RetainedMiniCollection )
        {
          // Couldn's retain this edge
          m_RetainedCost = itk::NumericTraits< float >::max();
        }
        else
        {
          m_RetainedCost = m_RetainedDataCost + m_RetainedAlphasCost + m_RetainedPositionCost;
        }



        // Calculate the cost of an edge collapse
        m_CollapsedCost = 0;
        m_CollapsedDataCost = 0;
        m_CollapsedAlphasCost = 0;
        m_CollapsedPositionCost = 0;
        std::set< AtlasMesh::CellIdentifier >  collapsedDisappearingCells;
        m_CollapsedMiniCollection =
          this->TryToCollapse( innerMiniCollection, outerMiniCollection, m_EdgeId,
                               m_CollapsedDataCost, m_CollapsedAlphasCost, m_CollapsedPositionCost,
                               collapsedDisappearingCells );
        if ( !m_CollapsedMiniCollection )
        {
          // Couldn't collapse this edge.
          m_CollapsedCost = itk::NumericTraits< float >::max();
        }
        else
        {
          m_CollapsedCost = m_CollapsedDataCost + m_CollapsedAlphasCost + m_CollapsedPositionCost;
        }


        // Evaluate which move is best
        std::cout << "                 retainedCost : " << m_RetainedCost
                  << "  (" << m_RetainedDataCost << " + " << m_RetainedAlphasCost
                  << " + " << m_RetainedPositionCost <<") " << std::endl;
        std::cout << "                 collapsedCost: " << m_CollapsedCost
                  << "  (" << m_CollapsedDataCost << " + " << m_CollapsedAlphasCost
                  << " + " << m_CollapsedPositionCost <<") " << std::endl;
        std::vector< float >  totalCosts;
        totalCosts.push_back( m_RetainedCost );
        totalCosts.push_back( m_CollapsedCost );
        float  minTotalCost = itk::NumericTraits< float >::max();
        int minTotalCostIndex = -1;
        for ( int i = 0; i < totalCosts.size(); i++ )
        {
          if ( totalCosts[ i ] < minTotalCost )
          {
            minTotalCost = totalCosts[ i ];
            minTotalCostIndex = i;
          }
        }

        if ( minTotalCostIndex == -1 )
        {
          std::cout << "Impossible configuration encountered at eget with id " << m_EdgeId << std::endl;
#if 1
          std::ostringstream  impossibleStream;
          impossibleStream << "impossible_" << m_EdgeId;
          m_Current->Write( impossibleStream.str().c_str() );
#endif
          //minTotalCostIndex = 0;
          continue;
        }



        // Do the best move
        AtlasMeshCollection::Pointer  newCurrent = AtlasMeshCollection::New();
        if ( minTotalCostIndex == 0 )
        {
          std::cout << "        => retaining edge is best solution" << std::endl;

          // this->ExpandRetain( const_cast< AtlasMeshCollection* >( m_Current.GetPointer() ), m_EdgeId,
          //                     m_RetainedMiniCollection,
          //                     newCurrent );


          // Look up the ids of the two points on the edge
          AtlasMesh::CellType::PointIdConstIterator  pointIt = m_Current->GetCells()->ElementAt( m_EdgeId )->PointIdsBegin();
          const AtlasMesh::PointIdentifier  edgePoint0Id = *pointIt;
          ++pointIt;
          const AtlasMesh::PointIdentifier  edgePoint1Id = *pointIt;

          // Remember the reference position of each of the two points
          newReferencePosition->InsertElement( edgePoint0Id,
                                               m_RetainedMiniCollection->GetReferencePosition()->ElementAt( edgePoint0Id ) );
          newReferencePosition->InsertElement( edgePoint1Id,
                                               m_RetainedMiniCollection->GetReferencePosition()->ElementAt( edgePoint1Id ) );

          // Remember the positions of each of the two points
          for ( unsigned int meshNumber = 0; meshNumber < m_Current->GetNumberOfMeshes(); meshNumber++ )
          {
            newPositions[ meshNumber ]->InsertElement( edgePoint0Id,
                m_RetainedMiniCollection->GetPositions()[ meshNumber ]->ElementAt( edgePoint0Id ) );
            newPositions[ meshNumber ]->InsertElement( edgePoint1Id,
                m_RetainedMiniCollection->GetPositions()[ meshNumber ]->ElementAt( edgePoint1Id ) );
          }

          // Remember the point parameters of each of the two points
          newPointParameters->InsertElement( edgePoint0Id,
                                             m_RetainedMiniCollection->GetPointParameters()->ElementAt( edgePoint0Id ) );
          newPointParameters->InsertElement( edgePoint1Id,
                                             m_RetainedMiniCollection->GetPointParameters()->ElementAt( edgePoint1Id ) );

        }
        else if ( minTotalCostIndex == 1 )
        {
          std::cout << "        => collapsing edge is best solution" << std::endl;

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

          // Look up the ids of the two points on the edge
          AtlasMesh::CellType::PointIdConstIterator  pointIt = m_Current->GetCells()->ElementAt( m_EdgeId )->PointIdsBegin();
          const AtlasMesh::PointIdentifier  edgePoint0Id = *pointIt;
          ++pointIt;
          const AtlasMesh::PointIdentifier  edgePoint1Id = *pointIt;


          // Insert them in the disappearingPointsLookupTable
          disappearingPointsLookupTable[ edgePoint1Id ] = edgePoint0Id;

          // Remember the reference position of the unified point
          newReferencePosition->InsertElement( edgePoint0Id,
                                               m_CollapsedMiniCollection->GetReferencePosition()->ElementAt( edgePoint0Id ) );

          // Remember the positions of the unified point
          for ( unsigned int meshNumber = 0; meshNumber < m_Current->GetNumberOfMeshes(); meshNumber++ )
          {
            newPositions[ meshNumber ]->InsertElement( edgePoint0Id,
                m_CollapsedMiniCollection->GetPositions()[ meshNumber ]->ElementAt( edgePoint0Id ) );
          }

          // Remember the point parameters of the unified point
          newPointParameters->InsertElement( edgePoint0Id,
                                             m_CollapsedMiniCollection->GetPointParameters()->ElementAt( edgePoint0Id ) );

          numberOfEdgeCollapses++;
        }




#if 0
        {
          std::ostringstream  outFileNameStream;
          outFileNameStream << "builderMeshCollection_Iteration" << m_IterationNumber << "_edge" << m_EdgeId << ".txt";
          m_Current->Write( outFileNameStream.str().c_str() );
        }
#endif

        // Next edge
        this->InvokeEvent( EndEdgeAnalysisEvent() );

      } // End loop over all edges


      std::cout << "Collapsed " << numberOfEdgeCollapses << " edges ("
                << static_cast< float >( numberOfEdgeCollapses ) /
                static_cast< float >( edgesToTry.size() ) * 100.0f
                << "% of all edges visited)" << std::endl;


      // Replace the current mesh with the new mesh
      m_Current = this->ApplyParallellMeshOperations( m_Current,
                  disappearingCells,
                  disappearingPointsLookupTable,
                  newReferencePosition,
                  newPositions,
                  newPointParameters ).GetPointer();


    } // End ICM iterations

    this->InvokeEvent( itk::EndEvent() );


  }

#elif 0


//
//
//
  void
  AtlasMeshBuilder
  ::Build( AtlasMeshCollection* explicitStartCollection )
  {

    if ( explicitStartCollection == 0 )
    {
      // Estimate with original mesh topology, using specified multi-resolution scheme
      AtlasParameterEstimator::LabelImageType::SizeType  labelImageSize =
        m_Estimator->GetLabelImage( 0 )->GetLargestPossibleRegion().GetSize();
      unsigned int  domainSize[ 3 ];
      domainSize[ 0 ] = static_cast< unsigned int >( labelImageSize[ 0 ] );
      domainSize[ 1 ] = static_cast< unsigned int >( labelImageSize[ 1 ] );
      domainSize[ 2 ] = static_cast< unsigned int >( labelImageSize[ 2 ] );

      for ( int upsamplingStepNumber = 0; upsamplingStepNumber <= m_NumberOfUpsamplingSteps; upsamplingStepNumber++ )
      {
        // Estimate
        m_Estimator->Estimate();

        // If this is not the final resolution yet, upsample mesh collection
        if ( upsamplingStepNumber != m_NumberOfUpsamplingSteps )
        {
          AtlasMeshCollection::Pointer  upsampledMeshCollection =
            m_Estimator->GetCurrentMeshCollection()->GetUpsampled();
          upsampledMeshCollection->SetK( m_InitialStiffnesses[ upsamplingStepNumber + 1 ] );
          m_Estimator->SetInitialMeshCollection( upsampledMeshCollection );

        }

      } // End loop over upsampling steps


    }
    else
    {
      // We have provided an explicit mesh already. Pass it on to the estimator
      m_Estimator->SetInitialMeshCollection( explicitStartCollection );
    }


    // Some initialization stuff
    m_Current = const_cast< AtlasMeshCollection* >( m_Estimator->GetCurrentMeshCollection() );

    this->InvokeEvent( itk::StartEvent() );



#if 0

    // ICM-like algorithm
    for ( m_IterationNumber = 0; m_IterationNumber < m_MaximumNumberOfIterations; m_IterationNumber++ )
    {
      this->InvokeEvent( itk::IterationEvent() );

      std::cout << "Iteration " << m_IterationNumber << std::endl;


      // Cache edges at this point
      std::vector< AtlasMesh::CellIdentifier >  edgesToTry;
      for ( AtlasMesh::CellsContainer::ConstIterator  cellIt = m_Current->GetCells()->Begin();
            cellIt != m_Current->GetCells()->End(); ++cellIt )
      {
        if ( cellIt.Value()->GetType() == AtlasMesh::CellType::LINE_CELL )
        {
          edgesToTry.push_back( cellIt.Index() );
        }
      }

      // Randomize the edges' order
      std::cout << "Randomizing edge order..." << std::endl;
      std::cout << "   number of edges before: " << edgesToTry.size() << std::endl;
      edgesToTry =  this->Permutate(  edgesToTry );
      std::cout << "   number of edges after: " << edgesToTry.size() << std::endl;

      // Loop over all edges
      unsigned int  progressCounter = 0;
      m_Progress = 0.0f;
      for ( std::vector< AtlasMesh::CellIdentifier >::const_iterator  edgeIt = edgesToTry.begin();
            edgeIt != edgesToTry.end(); ++edgeIt, ++progressCounter )
      {
        m_Progress = static_cast< float >( progressCounter ) / static_cast< float >( edgesToTry.size() );

        m_EdgeId = *edgeIt;

        this->InvokeEvent( StartEdgeAnalysisEvent() );

        this->AnalyzeEdge( *edgeIt );

        //std::ostringstream  outFileNameStream;
        //outFileNameStream << "builderMeshCollection_Iteration" << m_IterationNumber << "_edge" << m_EdgeId << ".txt";
        //m_Current->Write( outFileNameStream.str().c_str() );

        this->InvokeEvent( EndEdgeAnalysisEvent() );
      } // End loop over all edges

    } // End ICM iterations
#else

    itk::MultiThreader::Pointer  threader = itk::MultiThreader::New();
    threader->SetNumberOfThreads( 2 );
    int  numberOfThreads = threader->GetNumberOfThreads();
    std::cout << "Using " << numberOfThreads << " threads" << std::endl;


    // ICM-like algorithm
    for ( m_IterationNumber = 0; m_IterationNumber < m_MaximumNumberOfIterations; m_IterationNumber++ )
    {
      this->InvokeEvent( itk::IterationEvent() );

      std::cout << "Iteration " << m_IterationNumber << std::endl;

      int  numberOfEdgesAnalyzed = 0;
      while ( numberOfEdgesAnalyzed < 1000 )
      {
        // Get independent edges; one for each processor that's available
        ThreadStruct  str;
        str.m_Builder = this;
        //str.m_EdgesToTry = this->GetIndependentEdges( numberOfThreads );
        str.m_EdgesToTry = this->GetIndependentEdges();

        numberOfEdgesAnalyzed += str.m_EdgesToTry.size();
        // Print out what edges we're gonna analyze
        // std::cout << "str.m_EdgesToTry: [ ";
        // for ( std::vector< AtlasMesh::CellIdentifier >::const_iterator  it = str.m_EdgesToTry.begin();
        //       it != str.m_EdgesToTry.end(); ++it )
        //   {
        //   std::cout << *it << " ";
        //   }
        // std::cout << "]" << std::endl;


        // Do the job
        //threader->SetNumberOfThreads( str.m_EdgesToTry.size() );
        threader->SetSingleMethod( this->ThreaderCallback, &str );
        threader->SingleMethodExecute();
      }

      //
      // std::ostringstream  outFileNameStream;
      // outFileNameStream << "builderMeshCollection_ThreadedIteration" << m_IterationNumber << ".txt";
      // m_Current->Write( outFileNameStream.str().c_str() );

    } // End ICM iterations



#endif



    this->InvokeEvent( itk::EndEvent() );


  }


#elif 0


//
//
//
  void
  AtlasMeshBuilder
  ::Build( AtlasMeshCollection* explicitStartCollection )
  {

#if 0
    if ( explicitStartCollection == 0 )
    {
      // Estimate with original mesh topology, using specified multi-resolution scheme
      AtlasParameterEstimator::LabelImageType::SizeType  labelImageSize =
        m_Estimator->GetLabelImage( 0 )->GetLargestPossibleRegion().GetSize();
      unsigned int  domainSize[ 3 ];
      domainSize[ 0 ] = static_cast< unsigned int >( labelImageSize[ 0 ] );
      domainSize[ 1 ] = static_cast< unsigned int >( labelImageSize[ 1 ] );
      domainSize[ 2 ] = static_cast< unsigned int >( labelImageSize[ 2 ] );

      for ( int upsamplingStepNumber = 0; upsamplingStepNumber <= m_NumberOfUpsamplingSteps; upsamplingStepNumber++ )
      {
        // Estimate
        m_Estimator->Estimate();

        // If this is not the final resolution yet, upsample mesh collection
        if ( upsamplingStepNumber != m_NumberOfUpsamplingSteps )
        {
          AtlasMeshCollection::Pointer  upsampledMeshCollection =
            m_Estimator->GetCurrentMeshCollection()->GetUpsampled();
          upsampledMeshCollection->SetK( m_InitialStiffnesses[ upsamplingStepNumber + 1 ] );
          m_Estimator->SetInitialMeshCollection( upsampledMeshCollection );

        }

      } // End loop over upsampling steps

    }
    else
    {
      // We have provided an explicit mesh already. Pass it on to the estimator
      m_Estimator->SetInitialMeshCollection( explicitStartCollection );
    }


    // Some initialization stuff
    m_Current = const_cast< AtlasMeshCollection* >( m_Estimator->GetCurrentMeshCollection() );

#else
    if ( explicitStartCollection == 0 )
    {
      m_Mesher->Go();
      m_Current = const_cast< AtlasMeshCollection* >( m_Mesher->GetCurrentMeshCollection() );
    }
    else
    {
      m_Current = explicitStartCollection;
    }

#endif



    this->InvokeEvent( itk::StartEvent() );



#if 0

    // ICM-like algorithm
    for ( m_IterationNumber = 0; m_IterationNumber < m_MaximumNumberOfIterations; m_IterationNumber++ )
    {
      this->InvokeEvent( itk::IterationEvent() );

      std::cout << "Iteration " << m_IterationNumber << std::endl;


      // Get edges to analyze
      std::vector< AtlasMesh::CellIdentifier >  edgesToTry = this->GetIndependentEdges();

      // Some initialization stuff
      unsigned int  progressCounter = 0;
      m_Progress = 0.0f;
      std::map< AtlasMesh::PointIdentifier, AtlasMesh::PointIdentifier >   disappearingPointsLookupTable;
      AtlasMesh::PointsContainer::Pointer   newReferencePosition = AtlasMesh::PointsContainer::New();
      std::vector< AtlasMesh::PointsContainer::Pointer >   newPositions;
      int  numberOfEdgeCollapses = 0;
      for ( unsigned int meshNumber = 0; meshNumber < m_Current->GetNumberOfMeshes(); meshNumber++ )
      {
        newPositions.push_back( AtlasMesh::PointsContainer::New() );
      }
      AtlasMesh::PointDataContainer::Pointer  newPointParameters = AtlasMesh::PointDataContainer::New();
      std::set< AtlasMesh::CellIdentifier >   disappearingCells;


      // Loop over all edges
      for ( std::vector< AtlasMesh::CellIdentifier >::const_iterator  edgeIt = edgesToTry.begin();
            edgeIt != edgesToTry.end(); ++edgeIt, ++progressCounter )
      {
        m_Progress = static_cast< float >( progressCounter ) / static_cast< float >( edgesToTry.size() );

        m_EdgeId = *edgeIt;

        this->InvokeEvent( StartEdgeAnalysisEvent() );


        // Check if the edge is still present in the mesh.
        if ( !m_Current->GetCells()->IndexExists( m_EdgeId ) )
        {
          return;
        }


        // Get the mini collection surrounding the edge to try to collapse
        AtlasMeshCollection::ConstPointer  miniCollection = m_Current->GetRegionGrown( m_EdgeId, 1 ).GetPointer();

        // Retrieve the point ids on the edge
        AtlasMesh::CellType::PointIdConstIterator  pointIt = m_Current->GetCells()->ElementAt( m_EdgeId )->PointIdsBegin();
        const AtlasMesh::PointIdentifier  edgePoint0Id = *pointIt;
        ++pointIt;
        const AtlasMesh::PointIdentifier  edgePoint1Id = *pointIt;

        // Collect all the information about the collapsability of this edge
        this->AnalyzeEdgeFast( miniCollection, *edgeIt, edgePoint0Id, edgePoint1Id,
                               disappearingPointsLookupTable,
                               newReferencePosition, newPositions,
                               newPointParameters, disappearingCells );

        //std::ostringstream  outFileNameStream;
        //outFileNameStream << "builderMeshCollection_Iteration" << m_IterationNumber << "_edge" << m_EdgeId << ".txt";
        //m_Current->Write( outFileNameStream.str().c_str() );

        this->InvokeEvent( EndEdgeAnalysisEvent() );
      } // End loop over all edges


      // Apply the mesh surgery
      m_Current = this->ApplyParallellMeshOperations( m_Current,
                  disappearingCells,
                  disappearingPointsLookupTable,
                  newReferencePosition,
                  newPositions,
                  newPointParameters ).GetPointer();

    } // End ICM iterations
#else

    itk::MultiThreader::Pointer  threader = itk::MultiThreader::New();
    threader->SetNumberOfThreads( 2 );
    int  numberOfThreads = threader->GetNumberOfThreads();
    std::cout << "Using " << numberOfThreads << " threads" << std::endl;


    // ICM-like algorithm
    for ( m_IterationNumber = 0; m_IterationNumber < m_MaximumNumberOfIterations; m_IterationNumber++ )
    {
      this->InvokeEvent( itk::IterationEvent() );

      std::cout << "Iteration " << m_IterationNumber << std::endl;


#if 0
      int  numberOfEdgesAnalyzed = 0;
      while ( numberOfEdgesAnalyzed < 5000 )
      {
        // Get all the edges to loop over
        std::vector< AtlasMesh::CellIdentifier >  edgesToTry = this->GetIndependentEdges();
        //std::vector< AtlasMesh::CellIdentifier >  edgesToTry = this->GetIndependentEdges( 10 * numberOfThreads );
        numberOfEdgesAnalyzed += edgesToTry.size();
#else
      std::vector< AtlasMesh::CellIdentifier >  edges = this->GetRandomizedEdges();
      bool  neverRun = true;
      std::cout << "Inintial number of edges: " << edges.size() << std::endl;

      for ( std::vector< AtlasMesh::CellIdentifier >  edgesToTry =
              this->GetIndependentEdges( itk::NumericTraits< int >::max(), edges );
            ( edgesToTry.size() >= /* 100 */ 10 * numberOfThreads ) || neverRun;
            edgesToTry = this->GetIndependentEdges( itk::NumericTraits< int >::max(), edges ) )
      {

        if ( neverRun )
        {
          std::cout << "Going for first round with " << edgesToTry.size() << " edges to try" << std::endl;
        }
        else
        {
          std::cout << "Going for another round with " << edgesToTry.size() << " edges to try" << std::endl;
        }
        std::cout << "Remaining number of edges: " << edges.size() << std::endl;


        neverRun = false;
#endif

        // Build up the multi threading structure
        FastThreadStruct  str;
        str.m_Builder = this;
        for ( int threadNumber = 0; threadNumber < numberOfThreads; threadNumber++ )
        {
          // Construct an item for this thread
          //FastThreadStructItem*  item = new FastThreadStructItem;

          // Fill in the input stuff
          for ( int edgeNumber = threadNumber; edgeNumber < edgesToTry.size(); edgeNumber += numberOfThreads )
          {
            //
            const AtlasMesh::CellIdentifier  edgeId = edgesToTry[ edgeNumber ];

            // Get the mini collection surrounding the edge to try to collapse
            AtlasMeshCollection::ConstPointer  miniCollection = m_Current->GetRegionGrown( edgeId, 1 ).GetPointer();

            // Retrieve the point ids on the edge
            AtlasMesh::CellType::PointIdConstIterator  pointIt = m_Current->GetCells()->ElementAt( edgeId )->PointIdsBegin();
            const AtlasMesh::PointIdentifier  edgePoint0Id = *pointIt;
            ++pointIt;
            const AtlasMesh::PointIdentifier  edgePoint1Id = *pointIt;

            //
            str.m_Items[ threadNumber ].m_MiniCollections.push_back( miniCollection );
            str.m_Items[ threadNumber ].m_EdgeIds.push_back( edgeId );
            str.m_Items[ threadNumber ].m_EdgePoint0Ids.push_back( edgePoint0Id );
            str.m_Items[ threadNumber ].m_EdgePoint1Ids.push_back( edgePoint1Id );
          }

          // Initialize the output stuff
          str.m_Items[ threadNumber ].m_NewReferencePosition = AtlasMesh::PointsContainer::New();
          for ( unsigned int meshNumber = 0; meshNumber < m_Current->GetNumberOfMeshes(); meshNumber++ )
          {
            str.m_Items[ threadNumber ].m_NewPositions.push_back( AtlasMesh::PointsContainer::New() );
          }
          str.m_Items[ threadNumber ].m_NewPointParameters = AtlasMesh::PointDataContainer::New();

          // Now put the item in the multi threading structure
          //str->m_Items[ threadNumber ] = item;
        }

        // Do the job
        threader->SetSingleMethod( this->FastThreaderCallback, &str );
        threader->SingleMethodExecute();


        // Merge the output of each thread
        std::map< AtlasMesh::PointIdentifier, AtlasMesh::PointIdentifier >   disappearingPointsLookupTable;
        AtlasMesh::PointsContainer::Pointer   newReferencePosition = AtlasMesh::PointsContainer::New();
        std::vector< AtlasMesh::PointsContainer::Pointer >   newPositions;
        for ( unsigned int meshNumber = 0; meshNumber < m_Current->GetNumberOfMeshes(); meshNumber++ )
        {
          newPositions.push_back( AtlasMesh::PointsContainer::New() );
        }
        AtlasMesh::PointDataContainer::Pointer  newPointParameters = AtlasMesh::PointDataContainer::New();
        std::set< AtlasMesh::CellIdentifier >   disappearingCells;

        for ( int threadNumber = 0; threadNumber < numberOfThreads; threadNumber++ )
        {
          std::cout << "Merging output of thread number " << threadNumber << std::endl;

          // disappearingPointsLookupTable
          for ( std::map< AtlasMesh::PointIdentifier, AtlasMesh::PointIdentifier >::const_iterator  it =
                  str.m_Items[ threadNumber ].m_DisappearingPointsLookupTable.begin();
                it != str.m_Items[ threadNumber ].m_DisappearingPointsLookupTable.end(); ++it )
          {
            disappearingPointsLookupTable[ it->first ] = it->second;
          }

          // referencePosition
          for ( AtlasMesh::PointsContainer::ConstIterator  it = str.m_Items[ threadNumber ].m_NewReferencePosition->Begin();
                it != str.m_Items[ threadNumber ].m_NewReferencePosition->End(); ++it )
          {
            newReferencePosition->InsertElement( it.Index(), it.Value() );
          }

          // positions
          for ( unsigned int meshNumber = 0; meshNumber < m_Current->GetNumberOfMeshes(); meshNumber++ )
          {
            AtlasMesh::PointsContainer::ConstPointer  source =
              str.m_Items[ threadNumber ].m_NewPositions[ meshNumber ].GetPointer();
            AtlasMesh::PointsContainer::Pointer  target = newPositions[ meshNumber ];

            for ( AtlasMesh::PointsContainer::ConstIterator  it = source->Begin();
                  it != source->End(); ++it )
            {
              target->InsertElement( it.Index(), it.Value() );
            }

          }

          // point parameters
          for ( AtlasMesh::PointDataContainer::ConstIterator  it = str.m_Items[ threadNumber ].m_NewPointParameters->Begin();
                it != str.m_Items[ threadNumber ].m_NewPointParameters->End(); ++it )
          {
            newPointParameters->InsertElement( it.Index(), it.Value() );
          }

          // disappearing cells
          std::set< AtlasMesh::CellIdentifier >  totalDisappearingCells;
          std::set_union( disappearingCells.begin(), disappearingCells.end(),
                          str.m_Items[ threadNumber ].m_DisappearingCells.begin(),
                          str.m_Items[ threadNumber ].m_DisappearingCells.end(),
                          std::inserter( totalDisappearingCells, totalDisappearingCells.begin() ) );
          disappearingCells = totalDisappearingCells;

        }  // End loop over all thread outputs


        // Apply the mesh surgery
        m_Current = this->ApplyParallellMeshOperations( m_Current,
                    disappearingCells,
                    disappearingPointsLookupTable,
                    newReferencePosition,
                    newPositions,
                    newPointParameters ).GetPointer();

      }  // End test if we have already analyzed enough edges to start a new iteration



      //
      // std::ostringstream  outFileNameStream;
      // outFileNameStream << "builderMeshCollection_ThreadedIteration" << m_IterationNumber << ".txt";
      // m_Current->Write( outFileNameStream.str().c_str() );

    } // End ICM iterations



#endif



    this->InvokeEvent( itk::EndEvent() );

  }


#else


//
//
//
  void
  AtlasMeshBuilder
  ::Build( AtlasMeshCollection* explicitStartCollection )
  {

#if 0
    if ( explicitStartCollection == 0 )
    {
      // Estimate with original mesh topology, using specified multi-resolution scheme
      AtlasParameterEstimator::LabelImageType::SizeType  labelImageSize =
        m_Estimator->GetLabelImage( 0 )->GetLargestPossibleRegion().GetSize();
      unsigned int  domainSize[ 3 ];
      domainSize[ 0 ] = static_cast< unsigned int >( labelImageSize[ 0 ] );
      domainSize[ 1 ] = static_cast< unsigned int >( labelImageSize[ 1 ] );
      domainSize[ 2 ] = static_cast< unsigned int >( labelImageSize[ 2 ] );

      for ( int upsamplingStepNumber = 0; upsamplingStepNumber <= m_NumberOfUpsamplingSteps; upsamplingStepNumber++ )
      {
        // Estimate
        m_Estimator->Estimate();

        // If this is not the final resolution yet, upsample mesh collection
        if ( upsamplingStepNumber != m_NumberOfUpsamplingSteps )
        {
          AtlasMeshCollection::Pointer  upsampledMeshCollection =
            m_Estimator->GetCurrentMeshCollection()->GetUpsampled();
          upsampledMeshCollection->SetK( m_InitialStiffnesses[ upsamplingStepNumber + 1 ] );
          m_Estimator->SetInitialMeshCollection( upsampledMeshCollection );

        }

      } // End loop over upsampling steps


    }
    else
    {
      // We have provided an explicit mesh already. Pass it on to the estimator
      m_Estimator->SetInitialMeshCollection( explicitStartCollection );
    }

    // Some initialization stuff
    m_Current = const_cast< AtlasMeshCollection* >( m_Estimator->GetCurrentMeshCollection() );

#else

    if ( explicitStartCollection == 0 )
    {
      m_Mesher->Go();
      m_Current = const_cast< AtlasMeshCollection* >( m_Mesher->GetCurrentMeshCollection() );
    }
    else
    {
      m_Current = explicitStartCollection;
    }

#endif


    this->InvokeEvent( itk::StartEvent() );



    itk::MultiThreader::Pointer  threader = itk::MultiThreader::New();
#if 0
    if ( threader->GetNumberOfThreads() > 4 )
    {
      threader->SetNumberOfThreads( 4 );
    }
#endif
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


#endif




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
    std::cout << "Applying parallell mesh operations...." << std::endl;


    // Copy all positions of the points, except for the second point on the edge, which simply
    // disappears, and first point, which gets the previously determined new position
    std::cout << "Creating positions" << std::endl;
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
    std::cout << "Creating point parameters" << std::endl;

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
    std::cout << "Creating cells" << std::endl;
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


    std::cout << "...done!" << std::endl;

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
  ::TryToCollapse( const AtlasMeshCollection* innerMiniCollection, const AtlasMeshCollection* outerMiniCollection,
                   AtlasMesh::CellIdentifier  edgeId,
                   float& miniDataCost, float& miniAlphasCost, float& miniPositionCost,
                   std::set< AtlasMesh::CellIdentifier >&  disappearingCells ) const
  {
    //std::cout << "\n\n\n\n\n==================================================" << std::endl;
    //std::cout << "Trying to collapse edge " << edgeId << std::endl;



#if 0
    {
      std::ostringstream  meshStream;
      meshStream << "meshCollection_collapse" << edgeId << ".txt";
      meshCollection->Write( meshStream.str().c_str() );

      std::ostringstream  innerMiniStream;
      innerMiniStream << "innerMiniCollection_collapse" << edgeId << ".txt";
      innerMiniCollection->Write( innerMiniStream.str().c_str() );

      std::ostringstream  outerMiniStream;
      outerMiniStream << "outerMiniCollection_collapse" << edgeId << ".txt";
      outerMiniCollection->Write( outerMiniStream.str().c_str() );
    }
#endif

#if 0
    // Calculate the base line cost for the mini-mesh, so that the cost of the edge collapsed result
    // can be compared in order to calculate the gain
    float  dataCostOfOuterMiniCollection;
    float  alphasCostOfOuterMiniCollection;
    this->GetDataCostAndAlphasCost( outerMiniCollection, dataCostOfOuterMiniCollection, alphasCostOfOuterMiniCollection );
    const float positionCostOfOuterMiniCollection = this->GetPositionCost( outerMiniCollection );
#endif

    // Get the mesh with the edge collapsed.
    AtlasMesh::CellIdentifier  unifiedVertexId;
    AtlasMeshCollection::Pointer  outerChild;
    std::set< AtlasMesh::CellIdentifier >  disappearingCellsDummy;
    if ( !outerMiniCollection->GetCollapsed( edgeId, outerChild, disappearingCellsDummy, unifiedVertexId ) )
    {
      return 0;
    }
    AtlasMeshCollection::Pointer  innerChild;
    innerMiniCollection->GetCollapsed( edgeId, innerChild, disappearingCells, unifiedVertexId );


    const AtlasMesh::PointIdentifier  unifiedPointId =
      *( outerChild->GetCells()->ElementAt( unifiedVertexId )->PointIdsBegin() );

    // Set up the cost calculator
    kvl::AtlasMeshCollectionReferencePositionCost::ParametersType  initialPosition( 3 );
    initialPosition[ 0 ] = outerChild->GetReferencePosition()->ElementAt( unifiedPointId )[ 0 ];
    initialPosition[ 1 ] = outerChild->GetReferencePosition()->ElementAt( unifiedPointId )[ 1 ];
    initialPosition[ 2 ] = outerChild->GetReferencePosition()->ElementAt( unifiedPointId )[ 2 ];
    kvl::AtlasMeshCollectionReferencePositionCost::Pointer  costFunction =
      kvl::AtlasMeshCollectionReferencePositionCost::New();
    costFunction->SetLabelImages( m_LabelImages );
    costFunction->SetPointId( unifiedPointId );
    costFunction->SetInitialMeshCollections( innerChild, outerChild );
    kvl::AtlasMeshCollectionReferencePositionCost::ParametersType  optimalPosition = initialPosition;


    // Optimize the position of the unified vertex in the reference position
    const bool canMoveX = outerChild->GetPointParameters()->ElementAt( unifiedPointId ).m_CanMoveX;
    const bool canMoveY = outerChild->GetPointParameters()->ElementAt( unifiedPointId ).m_CanMoveY;
    const bool canMoveZ = outerChild->GetPointParameters()->ElementAt( unifiedPointId ).m_CanMoveZ;
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
      //optimizer->SetAbsolutePrecisionBrent( m_PowellAbsolutePrecision );
      optimizer->SetStepTolerance( m_PowellAbsolutePrecision );
      //optimizer->SetAbsolutePrecision( m_PowellAbsolutePrecision );
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


#if 0
    if ( !canMoveX && !canMoveY )
    {
      std::cout << "Having a special case here. We didn't optimize, and the result is:" << std::endl;
      std::cout << "       miniDataCost     : " << miniDataCost << std::endl;
      std::cout << "       miniAlphasCost   : " << miniAlphasCost << std::endl;
      std::cout << "       miniPositionCost : " << miniPositionCost << std::endl;
      innerMiniCollection->Write( "specialCase_innerMiniCollection.txt" );
      outerMiniCollection->Write( "specialCase_outerMiniCollection.txt" );
      innerChild->Write( "specialCase_innerChild.txt" );
      outerChild->Write( "specialCase_outerChild.txt" );
      costFunction->GetCostCalculationMeshCollection()->Write( "specialCaseHere.txt" );
      exit( -1 );
    }

#endif



#if 0
    std::ostringstream  childStream;
    childStream << "debugChildCollection_iteration" <<  m_EdgeCollapseNumber
                << "_edge" << edgeId;
    child->Write( childStream.str().c_str() );

    std::ostringstream  dataAndAlphasCostStream;
    dataAndAlphasCostStream << "debugDataAndAlphasCost_iteration" <<  m_EdgeCollapseNumber
                            << "_edge" << edgeId;

    std::ofstream  out( dataAndAlphasCostStream.str().c_str(), std::ios::app );
    if ( !out.bad() )
    {
      out << "miniCollection: " << std::endl;
      out << "                   dataCost: " << dataCostOfMiniCollection << std::endl;
      out << "                 alphasCost: " << alphasCostOfMiniCollection << std::endl;
      out << "               positionCost: " << positionCostOfMiniCollection << std::endl;
      out << std::endl;
      out << "child: " << std::endl;
      out << "                   dataCost: " << childDataCost << std::endl;
      out << "                 alphasCost: " << childAlphasCost << std::endl;
      out << "               positionCost: " << childPositionCost << std::endl;
      out << "---------------------------------------------" << std::endl;
      out << "gain: " << std::endl;
      out << "                   dataCost: " << dataCostOfMiniCollection - childDataCost << std::endl;
      out << "                 alphasCost: " << alphasCostOfMiniCollection - childAlphasCost << std::endl;
      out << "               positionCost: " << positionCostOfMiniCollection - childPositionCost << std::endl;
    }
#endif


#if 0
#if 1
    dataGain =  dataCostOfOuterMiniCollection - outerChildDataCost;
    alphasGain = alphasCostOfOuterMiniCollection - outerChildAlphasCost;
    positionGain = positionCostOfOuterMiniCollection - outerChildPositionCost;

#else
    // Relative gains!
    const float totalCost =
      ( dataCostOfOuterMiniCollection + alphasCostOfOuterMiniCollection + positionCostOfOuterMiniCollection +
        outerChildDataCost + outerChildAlphasCost + outerChildPositionCost ) / 2;
    dataGain =  ( dataCostOfOuterMiniCollection - outerChildDataCost ) / totalCost;
    alphasGain = ( alphasCostOfOuterMiniCollection - outerChildAlphasCost ) / totalCost;
    positionGain = ( positionCostOfOuterMiniCollection - outerChildPositionCost ) / totalCost ;


#endif
#endif


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
    costFunction->SetPointId( unifiedPointId );
    costFunction->SetInitialMeshCollection( child );
    kvl::AtlasMeshCollectionFastReferencePositionCost::ParametersType  optimalPosition = initialPosition;


#if 1

    if ( !optimizeUnifiedPointReferencePosition )
    {
      std::cout << "NOT optimizing unified reference position" << std::endl;
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
    std::cout << "\n\n\n\n\n==================================================" << std::endl;
    std::cout << "Trying to expand edge collapse for edge " << edgeId << std::endl;




    // Get the mesh with the edge collapsed.
    if ( !meshCollection->GetCollapsed( edgeId, result, *disappearingCells, unifiedVertexId ) )
    {
      itkExceptionMacro( "Couldn't initialize positions with those estimated in mini-mesh!" );
    }


    // Copy the positions estimated in the mini mesh to the collapsed collection
    std::vector< AtlasMesh::PointsContainer::Pointer >  childPositions =  optimizedOuterChild->GetPositions();
    std::cout << "Copying positions from mini mesh collection" << std::endl;
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

#if 0
    std::cout << "\n\n\n\n!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << std::endl;
    this->DebugOn();
    const float  positionCostOfOuterMiniCollection = this->GetPositionCost( outerMiniCollection );
    const float  positionCostOfOptimizedOuterChild = this->GetPositionCost( optimizedOuterChild );
    const float  positionCostOfMeshCollection = this->GetPositionCost( meshCollection );
    const float  positionCostOfResult = this->GetPositionCost( result );
    std::cout << "Position cost difference on mini mesh: " << positionCostOfOuterMiniCollection
              << " - " <<  positionCostOfOptimizedOuterChild
              << " = " << positionCostOfOuterMiniCollection - positionCostOfOptimizedOuterChild << std::endl;
    std::cout << "Position cost difference on full mesh: " << positionCostOfMeshCollection
              << " - " <<  positionCostOfResult
              << " = " << positionCostOfMeshCollection - positionCostOfResult << std::endl;
    std::cout << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n\n\n\n" << std::endl;

    outerMiniCollection->Write( "outerMiniCollection.txt" );
    optimizedOuterChild->Write( "optimizedOuterChild.txt" );
    meshCollection->Write( "meshCollection.txt" );
    result->Write( "result.txt" );
    exit( 0 );
#endif



  }



#if 0
//
//
//
  AtlasMeshCollection::Pointer
  AtlasMeshBuilder
  ::TryToSplit( const AtlasMeshCollection* meshCollection, AtlasMesh::CellIdentifier  edgeId,
                float& miniDataCost, float& miniAlphasCost, float& miniPositionCost,
                AtlasMeshCollection::Pointer& result ) const
  {
    //std::cout << "\n\n\n\n\n==================================================" << std::endl;
    //std::cout << "Trying to split edge " << edgeId << std::endl;


    // Get the mini collection surrounding the edge to try to split
    AtlasMeshCollection::ConstPointer  outerMiniCollection = 0;
    AtlasMeshCollection::ConstPointer  innerMiniCollection = 0;
    if ( !m_AllowMiniMeshSpeedUp )
    {
      // The mini-mesh is the entire mesh
      outerMiniCollection = meshCollection;
      innerMiniCollection = meshCollection;
    }
    else
    {
      // The mini-mesh is an area surrounding the edge to collapse
      innerMiniCollection = meshCollection->GetRegionGrown( edgeId, m_MiniMeshRadius );
      outerMiniCollection = meshCollection->GetRegionGrown( edgeId, m_MiniMeshRadius+1 );
    }

#if 0
    {
      std::ostringstream  meshStream;
      meshStream << "meshCollection_split" << edgeId << ".txt";
      meshCollection->Write( meshStream.str().c_str() );

      std::ostringstream  innerMiniStream;
      innerMiniStream << "innerMiniCollection_split" << edgeId << ".txt";
      innerMiniCollection->Write( innerMiniStream.str().c_str() );

      std::ostringstream  outerMiniStream;
      outerMiniStream << "outerMiniCollection_split" << edgeId << ".txt";
      outerMiniCollection->Write( outerMiniStream.str().c_str() );
    }
#endif

#if 0
    {
      m_Estimator->SetInitialMeshCollection( const_cast< AtlasMeshCollection* >( outerMiniCollection.GetPointer() ) );
      m_Estimator->Estimate();
      outerMiniCollection = m_Estimator->GetCurrentMeshCollection();
    }
#endif

#if 0
    // Calculate the base line cost for the mini-mesh, so that the cost of the edge split result
    // can be compared in order to calculate the gain
    float  dataCostOfOuterMiniCollection;
    float  alphasCostOfOuterMiniCollection;
    this->GetDataCostAndAlphasCost( outerMiniCollection, dataCostOfOuterMiniCollection, alphasCostOfOuterMiniCollection );
    const float positionCostOfOuterMiniCollection = this->GetPositionCost( outerMiniCollection );
#endif

    // Get the mesh with the edge split.
    AtlasMesh::CellsContainer::ConstIterator  lastCellIt = meshCollection->GetCells()->End();
    lastCellIt--;
    AtlasMesh::PointsContainer::ConstIterator  lastPointIt = meshCollection->GetReferencePosition()->End();
    lastPointIt--;
    const AtlasMesh::CellIdentifier  newVertexId = lastCellIt.Index() + 1;
    const AtlasMesh::PointIdentifier  newPointId = lastPointIt.Index() + 1;
    //std::cout << "Splitting edge " << edgeId << " (newVertexId: " << newVertexId
    //          << ", newPointId: " << newPointId<< ")" << std::endl;
    AtlasMeshCollection::Pointer  outerChild = outerMiniCollection->GetEdgeSplitted( edgeId, newVertexId, newPointId );
    if ( !outerChild )
    {
      return 0;
    }
    AtlasMeshCollection::Pointer  innerChild = innerMiniCollection->GetEdgeSplitted( edgeId, newVertexId, newPointId );

#if 0
    {
      std::ostringstream  innerChildStream;
      innerChildStream << "innerChild_" << edgeId << ".txt";
      innerChild->Write( innerChildStream.str().c_str() );

      std::ostringstream  outerChildStream;
      outerChildStream << "outerChild_" << edgeId << ".txt";
      outerChild->Write( outerChildStream.str().c_str() );

      /*  char  dummy;
        std::cin >> dummy;*/
    }
#endif

    // Set up the cost calculator
    kvl::AtlasMeshCollectionReferencePositionCost::ParametersType  initialPosition( 2 );
    initialPosition[ 0 ] = outerChild->GetReferencePosition()->ElementAt( newPointId )[ 0 ];
    initialPosition[ 1 ] = outerChild->GetReferencePosition()->ElementAt( newPointId )[ 1 ];
    kvl::AtlasMeshCollectionReferencePositionCost::Pointer  costFunction =
      kvl::AtlasMeshCollectionReferencePositionCost::New();
    costFunction->SetLabelImages( m_LabelImages );
    costFunction->SetPointId( newPointId );
    costFunction->SetInitialMeshCollections( innerChild, outerChild );
    kvl::AtlasMeshCollectionReferencePositionCost::ParametersType  optimalPosition = initialPosition;


    // Optimize the position of the unified vertex in the reference position
    const bool canMoveX = outerChild->GetPointParameters()->ElementAt( newPointId ).m_CanMoveX;
    const bool canMoveY = outerChild->GetPointParameters()->ElementAt( newPointId ).m_CanMoveY;
    if ( canMoveX || canMoveY )
    {

      // Decide what to optimize
      itk::PowellOptimizer::ParameterOrderType  parameterOrder( 2 );
      if ( canMoveX && !canMoveY )
      {
        parameterOrder[ 0 ] = 1;
        parameterOrder[ 1 ] = 0;
      }
      else if ( !canMoveX && canMoveY )
      {
        parameterOrder[ 0 ] = 0;
        parameterOrder[ 1 ] = 1;
      }
      else
      {
        parameterOrder[ 0 ] = 1;
        parameterOrder[ 1 ] = 2;
      }


      // Optimize the reference position
      //std::cout << "Optimizing the position of the unified point " << unifiedPointId << std::endl;
      itk::PowellOptimizer::Pointer  optimizer = itk::PowellOptimizer::New();
      optimizer->SetCostFunction( costFunction );
      optimizer->SetInitialPosition( initialPosition );
      optimizer->SetAbsolutePrecisionBrent( m_PowellAbsolutePrecision );
      optimizer->SetAbsolutePrecision( m_PowellAbsolutePrecision );
      optimizer->SetParameterOrder( parameterOrder );
      optimizer->StartOptimization();

      // Retrieve the optimal reference position
      if ( optimizer->GetCurrentValue() != itk::NumericTraits< float >::max() )
      {
        optimalPosition = optimizer->GetCurrentPosition();
        //std::cout << "Changed reference position for the new point " << newPointId
        //          << " from " << initialPosition << " to " << optimalPosition << std::endl;
      }
      else
      {
        return 0;
      }



    }  // End test if point can move



    // Get the cost at the optimal position
    if ( !costFunction->GetValue( optimalPosition, miniDataCost, miniAlphasCost, miniPositionCost ) )
    {
      return 0;
    }

#if 0
#if 1
    dataGain =  dataCostOfOuterMiniCollection - outerChildDataCost;
    alphasGain = alphasCostOfOuterMiniCollection - outerChildAlphasCost;
    positionGain = positionCostOfOuterMiniCollection - outerChildPositionCost;

#else
    // Relative gains!
    const float totalCost =
      ( dataCostOfOuterMiniCollection + alphasCostOfOuterMiniCollection + positionCostOfOuterMiniCollection +
        outerChildDataCost + outerChildAlphasCost + outerChildPositionCost ) / 2;
    dataGain =  ( dataCostOfOuterMiniCollection - outerChildDataCost ) / totalCost;
    alphasGain = ( alphasCostOfOuterMiniCollection - outerChildAlphasCost ) / totalCost;
    positionGain = ( positionCostOfOuterMiniCollection - outerChildPositionCost ) / totalCost ;

#endif
#endif

    if ( result )
    {
      // Get the mesh with the edge collapsed.
      result = meshCollection->GetEdgeSplitted( edgeId, newVertexId, newPointId );
      if ( !result )
      {
        itkExceptionMacro( "Couldn't initialize positions with those estimated in mini-mesh!" );
      }

      // Get the optimized outer child mesh from the costFunction
      AtlasMeshCollection::ConstPointer  optimizedOuterChild = costFunction->GetCostCalculationMeshCollection();

#if 0
      {
        std::ostringstream  outerMiniStream;
        outerMiniStream << "outerMiniCollection_split" << edgeId << ".txt";
        outerMiniCollection->Write( outerMiniStream.str().c_str() );

        std::ostringstream  optimizedOuterChildStream;
        optimizedOuterChildStream << "optimizedOuterChild_split" << edgeId << ".txt";
        optimizedOuterChild->Write( optimizedOuterChildStream.str().c_str() );
      }
#endif


      // Copy the positions estimated in the mini mesh to the splitted collection
      std::vector< AtlasMesh::PointsContainer::Pointer >  childPositions =  optimizedOuterChild->GetPositions();
      //std::cout << "Copying positions from mini mesh collection" << std::endl;
      for ( unsigned int meshNumber = 0; meshNumber < childPositions.size(); meshNumber++ )
      {
        AtlasMesh::PointsContainer::ConstIterator  positionIt = childPositions[ meshNumber ]->Begin();
        while ( positionIt !=  childPositions[ meshNumber ]->End() )
        {
          result->GetPositions()[ meshNumber ]->ElementAt( positionIt.Index() ) =  positionIt.Value();

          ++positionIt;
        }

      }



      // Also copy the reference position from the child to the original mesh after splitting the edge
      //std::cout << "Copying the reference position from the mini-mesh to the result mesh" << std::endl;
      result->GetReferencePosition()->ElementAt( newPointId ) =
        optimizedOuterChild->GetReferencePosition()->ElementAt( newPointId );

      // Force the triangle area parameters to be updated!!!
      result->SetK( result->GetK() );


      // First copy the alphas from the original mesh, and then overwrite the alphas within the child
      for ( AtlasMesh::PointDataContainer::ConstIterator  paramIt = meshCollection->GetPointParameters()->Begin();
            paramIt != meshCollection->GetPointParameters()->End();
            ++paramIt )
      {
        result->GetPointParameters()->ElementAt( paramIt.Index() ).m_Alphas = paramIt.Value().m_Alphas;
      }
      for ( AtlasMesh::PointDataContainer::ConstIterator  paramIt = optimizedOuterChild->GetPointParameters()->Begin();
            paramIt != optimizedOuterChild->GetPointParameters()->End();
            ++paramIt )
      {
        result->GetPointParameters()->ElementAt( paramIt.Index() ).m_Alphas = paramIt.Value().m_Alphas;
      }


    }

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
  AtlasMeshCollection::Pointer
  AtlasMeshBuilder
  ::TryToSwap( const AtlasMeshCollection* meshCollection, AtlasMesh::CellIdentifier  edgeId,
               float& miniDataCost, float& miniAlphasCost, float& miniPositionCost,
               AtlasMeshCollection::Pointer& result ) const
  {
    //std::cout << "\n\n\n\n\n==================================================" << std::endl;
    //std::cout << "Trying to swap edge " << edgeId << std::endl;


    // Get the mini collection surrounding the edge to try to swap
    AtlasMeshCollection::ConstPointer  outerMiniCollection = 0;
    AtlasMeshCollection::ConstPointer  innerMiniCollection = 0;
    if ( !m_AllowMiniMeshSpeedUp )
    {
      // The mini-mesh is the entire mesh
      outerMiniCollection = meshCollection;
      innerMiniCollection = meshCollection;
    }
    else
    {
      // The mini-mesh is an area surrounding the edge to collapse
      innerMiniCollection = meshCollection->GetRegionGrown( edgeId, m_MiniMeshRadius );
      outerMiniCollection = meshCollection->GetRegionGrown( edgeId, m_MiniMeshRadius+1 );
    }

#if 0
    {
      std::ostringstream  meshStream;
      meshStream << "meshCollection_swap" << edgeId << ".txt";
      meshCollection->Write( meshStream.str().c_str() );

      //std::ostringstream  innerMiniStream;
      //innerMiniStream << "innerMiniCollection_swap" << edgeId << ".txt";
      //innerMiniCollection->Write( innerMiniStream.str().c_str() );

      std::ostringstream  outerMiniStream;
      outerMiniStream << "outerMiniCollection_swap" << edgeId << ".txt";
      outerMiniCollection->Write( outerMiniStream.str().c_str() );
    }
#endif


#if 0
    // Calculate the base line cost for the mini-mesh, so that the cost of the edge split result
    // can be compared in order to calculate the gain
    float  dataCostOfOuterMiniCollection;
    float  alphasCostOfOuterMiniCollection;
    this->GetDataCostAndAlphasCost( outerMiniCollection, dataCostOfOuterMiniCollection, alphasCostOfOuterMiniCollection );
    const float positionCostOfOuterMiniCollection = this->GetPositionCost( outerMiniCollection );
#endif


    // Get the mesh with the edge swapped.
    AtlasMesh::CellsContainer::ConstIterator  lastCellIt = meshCollection->GetCells()->End();
    lastCellIt--;
    AtlasMesh::PointsContainer::ConstIterator  lastPointIt = meshCollection->GetReferencePosition()->End();
    lastPointIt--;
    const AtlasMesh::CellIdentifier  newVertexId = lastCellIt.Index() + 1;
    const AtlasMesh::PointIdentifier  newPointId = lastPointIt.Index() + 1;
    //std::cout << "Swapping edge " << edgeId << " (newVertexId: " << newVertexId
    //          << ", newPointId: " << newPointId<< ")" << std::endl;
    AtlasMesh::CellIdentifier  newEdgeId;
    AtlasMeshCollection::Pointer  outerChild =
      outerMiniCollection->GetEdgeSwapped( edgeId, newVertexId, newPointId, newEdgeId );
    if ( !outerChild )
    {
      return 0;
    }
    //AtlasMeshCollection::Pointer  innerChild =
    //          innerMiniCollection->GetEdgeSwapped( edgeId, newVertexId, newPointId, newEdgeId );

#if 0
    {
      //std::ostringstream  innerChildStream;
      //innerChildStream << "innerChild_swap" << edgeId << ".txt";
      //innerChild->Write( innerChildStream.str().c_str() );

      std::ostringstream  outerChildStream;
      outerChildStream << "outerChild_swap" << edgeId << ".txt";
      outerChild->Write( outerChildStream.str().c_str() );

      /*  char  dummy;
        std::cin >> dummy;*/
    }
#endif



    // Retrieve the point ids of the vertices of the newly created edge
    AtlasMesh::CellType::PointIdConstIterator  pointIt =
      outerChild->GetCells()->ElementAt( newEdgeId )->PointIdsBegin();
    const AtlasMesh::PointIdentifier  point0Id = *pointIt;
    ++pointIt;
    const AtlasMesh::PointIdentifier  point1Id = *pointIt;


    // Retrieve the vertex ids of the vertices of this edge, by looping over all vertex cells and
    // checking their point id
    AtlasMesh::CellIdentifier  vertex0Id;
    AtlasMesh::CellIdentifier  vertex1Id;
    for ( AtlasMesh::CellsContainer::ConstIterator  cellIt = outerChild->GetCells()->Begin();
          cellIt != outerChild->GetCells()->End(); ++cellIt )
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
    if ( !this->OptimizeReferencePosition( outerChild, vertex0Id ) )
    {
      return 0;
    }
    if ( !this->OptimizeReferencePosition( outerChild, vertex1Id ) )
    {
      return 0;
    }

#if 0
    {
      //std::ostringstream  innerChildStream;
      //innerChildStream << "innerChild_swap" << edgeId << ".txt";
      //innerChild->Write( innerChildStream.str().c_str() );

      std::ostringstream  outerChildStream;
      outerChildStream << "outerChild_afterOptimization_swap" << edgeId << ".txt";
      outerChild->Write( outerChildStream.str().c_str() );

      char  dummy;
      std::cin >> dummy;
    }
#endif

    // Measure the cost after optimization
    this->GetDataCostAndAlphasCost( outerChild, miniDataCost, miniAlphasCost );
    miniPositionCost = this->GetPositionCost( outerChild );


#if 0
    // Set up the cost calculator; we'll abuse it here to simply calculate for us the optimized outer child
    // and its cost
    const AtlasMesh::PointIdentifier  pretentionPointId = innerChild->GetReferencePosition()->Begin().Index();
    kvl::AtlasMeshCollectionReferencePositionCost::ParametersType  pretentionPosition( 2 );
    pretentionPosition[ 0 ] = innerChild->GetReferencePosition()->Begin().Value()[ 0 ];
    pretentionPosition[ 1 ] = innerChild->GetReferencePosition()->Begin().Value()[ 1 ];
    kvl::AtlasMeshCollectionReferencePositionCost::Pointer  costFunction =
      kvl::AtlasMeshCollectionReferencePositionCost::New();
    costFunction->SetLabelImages( m_LabelImages );
    costFunction->SetPointId( pretentionPointId );
    costFunction->SetInitialMeshCollections( innerChild, outerChild );


    // Get the cost at the optimal position
    if ( !costFunction->GetValue( pretentionPosition, miniDataCost, miniAlphasCost, miniPositionCost ) )
    {
      return false;
    }

#if 0
#if 1
    dataGain =  dataCostOfOuterMiniCollection - outerChildDataCost;
    alphasGain = alphasCostOfOuterMiniCollection - outerChildAlphasCost;
    positionGain = positionCostOfOuterMiniCollection - outerChildPositionCost;

#else
    // Relative gains!
    const float totalCost =
      ( dataCostOfOuterMiniCollection + alphasCostOfOuterMiniCollection + positionCostOfOuterMiniCollection +
        outerChildDataCost + outerChildAlphasCost + outerChildPositionCost ) / 2;
    dataGain =  ( dataCostOfOuterMiniCollection - outerChildDataCost ) / totalCost;
    alphasGain = ( alphasCostOfOuterMiniCollection - outerChildAlphasCost ) / totalCost;
    positionGain = ( positionCostOfOuterMiniCollection - outerChildPositionCost ) / totalCost ;

#endif
#endif

#endif

    if ( result )
    {
      // Get the mesh with the edge swapped.
      result = meshCollection->GetEdgeSwapped( edgeId, newVertexId, newPointId, newEdgeId );
      if ( !result )
      {
        itkExceptionMacro( "Couldn't initialize positions with those estimated in mini-mesh!" );
      }


      // Copy the positions estimated in the mini mesh to the swapped collection
      std::vector< AtlasMesh::PointsContainer::Pointer >  childPositions =  outerChild->GetPositions();
      //std::cout << "Copying positions from mini mesh collection" << std::endl;
      for ( unsigned int meshNumber = 0; meshNumber < childPositions.size(); meshNumber++ )
      {
        AtlasMesh::PointsContainer::ConstIterator  positionIt = childPositions[ meshNumber ]->Begin();
        while ( positionIt !=  childPositions[ meshNumber ]->End() )
        {
          result->GetPositions()[ meshNumber ]->ElementAt( positionIt.Index() ) =  positionIt.Value();

          ++positionIt;
        }

      }

      // Also copy the reference position from the mini mesh to the original mesh
      //std::cout << "Copying the reference positions from the mini-mesh to the result mesh" << std::endl;
      result->GetReferencePosition()->ElementAt( point0Id ) =
        outerChild->GetReferencePosition()->ElementAt( point0Id );
      result->GetReferencePosition()->ElementAt( point1Id ) =
        outerChild->GetReferencePosition()->ElementAt( point1Id );

      // Force the triangle area parameters to be updated!!!
      result->SetK( meshCollection->GetK() );

      // First copy the alphas from the original mesh, and then overwrite the alphas within the child
      for ( AtlasMesh::PointDataContainer::ConstIterator  paramIt = meshCollection->GetPointParameters()->Begin();
            paramIt != meshCollection->GetPointParameters()->End();
            ++paramIt )
      {
        result->GetPointParameters()->ElementAt( paramIt.Index() ).m_Alphas = paramIt.Value().m_Alphas;
      }
      for ( AtlasMesh::PointDataContainer::ConstIterator  paramIt = outerChild->GetPointParameters()->Begin();
            paramIt != outerChild->GetPointParameters()->End();
            ++paramIt )
      {
        result->GetPointParameters()->ElementAt( paramIt.Index() ).m_Alphas = paramIt.Value().m_Alphas;
      }


    }

    // Return
    if ( std::isnan( miniDataCost + miniAlphasCost + miniPositionCost ) )
    {
      return 0;
    }
    else
    {
      return outerChild;
    }

  }
#endif






//
//
//
  AtlasMeshCollection::Pointer
  AtlasMeshBuilder
  ::TryToRetain( const AtlasMeshCollection* innerMiniCollectionConst, const AtlasMeshCollection* outerMiniCollectionConst,
                 AtlasMesh::CellIdentifier  edgeId,
                 float& miniDataCost, float& miniAlphasCost, float& miniPositionCost )
  {

    //std::cout << "\n\n\n\n\n==================================================" << std::endl;
    //std::cout << "Trying to retain edge " << edgeId << std::endl;

#if 0
    // Calculate the base line cost for the mini-mesh, so that the cost of the edge retain result
    // can be compared in order to calculate the gain
    float  dataCostOfOuterMiniCollection;
    float  alphasCostOfOuterMiniCollection;
    this->GetDataCostAndAlphasCost( outerMiniCollection, dataCostOfOuterMiniCollection, alphasCostOfOuterMiniCollection );
    const float positionCostOfOuterMiniCollection = this->GetPositionCost( outerMiniCollection );
#endif


    // Make a non-const copy of mesh collections.
#if 0
// TODO: do this differently
    innerMiniCollectionConst->Write( "inner.txt" );
    outerMiniCollectionConst->Write( "outer.txt" );

    AtlasMeshCollection::Pointer  innerMiniCollection = AtlasMeshCollection::New();
    innerMiniCollection->Read( "inner.txt" );
    AtlasMeshCollection::Pointer  outerMiniCollection = AtlasMeshCollection::New();
    outerMiniCollection->Read( "outer.txt" );
#else
    AtlasMeshCollection::Pointer  innerMiniCollection = this->GetFakeCopy( innerMiniCollectionConst );
    AtlasMeshCollection::Pointer  outerMiniCollection = this->GetFakeCopy( outerMiniCollectionConst );

    // innerMiniCollection->Write( "inner.txt" );
    // outerMiniCollection->Write( "outer.txt" );
#endif


    // Retrieve the point ids of the vertices of this edge
    AtlasMesh::CellType::PointIdConstIterator  pointIt =
      innerMiniCollection->GetCells()->ElementAt( edgeId )->PointIdsBegin();
    const AtlasMesh::PointIdentifier  point0Id = *pointIt;
    ++pointIt;
    const AtlasMesh::PointIdentifier  point1Id = *pointIt;


    // Retrieve the vertex ids of the vertices of this edge, by looping over all vertex cells and
    // checking their point id
    AtlasMesh::CellIdentifier  vertex0Id;
    AtlasMesh::CellIdentifier  vertex1Id;
    for ( AtlasMesh::CellsContainer::ConstIterator  cellIt = innerMiniCollection->GetCells()->Begin();
          cellIt != innerMiniCollection->GetCells()->End(); ++cellIt )
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
    if ( !this->OptimizeReferencePosition( innerMiniCollection, vertex0Id, false ) )
    {
      return 0;
    }
    if ( !this->OptimizeReferencePosition( innerMiniCollection, vertex1Id, false ) )
    {
      return 0;
    }




    // Expand the inner mini mesh to the outer mini mesh

    // Copy the positions estimated in the inner mini mesh to the outer mini mesh
    std::vector< AtlasMesh::PointsContainer::Pointer >  innerPositions =  innerMiniCollection->GetPositions();
    for ( unsigned int meshNumber = 0; meshNumber < innerPositions.size(); meshNumber++ )
    {
      for ( AtlasMesh::PointsContainer::ConstIterator  positionIt = innerPositions[ meshNumber ]->Begin();
            positionIt !=  innerPositions[ meshNumber ]->End(); ++positionIt )
      {
        outerMiniCollection->GetPositions()[ meshNumber ]->ElementAt( positionIt.Index() ) =  positionIt.Value();
      }

    }

    // Also copy the reference position from the mini mesh to the original mesh
    outerMiniCollection->GetReferencePosition()->ElementAt( point0Id ) =
      innerMiniCollection->GetReferencePosition()->ElementAt( point0Id );
    outerMiniCollection->GetReferencePosition()->ElementAt( point1Id ) =
      innerMiniCollection->GetReferencePosition()->ElementAt( point1Id );

    // Force the triangle area parameters to be updated!!!
    outerMiniCollection->SetK( outerMiniCollection->GetK() );


    // Copy the alphas from the mini mesh to the whole mesh
    for ( AtlasMesh::PointDataContainer::ConstIterator  paramIt = innerMiniCollection->GetPointParameters()->Begin();
          paramIt != innerMiniCollection->GetPointParameters()->End();
          ++paramIt )
    {
      outerMiniCollection->GetPointParameters()->ElementAt( paramIt.Index() ).m_Alphas = paramIt.Value().m_Alphas;
    }





    // Measure the cost after optimization
    this->GetDataCostAndAlphasCost( outerMiniCollection, miniDataCost, miniAlphasCost );
    miniPositionCost = this->GetPositionCost( outerMiniCollection );

#if 0
#if 1
    dataGain =  dataCostOfOuterMiniCollection - outerChildDataCost;
    alphasGain = alphasCostOfOuterMiniCollection - outerChildAlphasCost;
    positionGain = positionCostOfOuterMiniCollection - outerChildPositionCost;

#else
    // Relative gains!
    const float totalCost =
      ( dataCostOfOuterMiniCollection + alphasCostOfOuterMiniCollection + positionCostOfOuterMiniCollection +
        outerChildDataCost + outerChildAlphasCost + outerChildPositionCost ) / 2;
    dataGain =  ( dataCostOfOuterMiniCollection - outerChildDataCost ) / totalCost;
    alphasGain = ( alphasCostOfOuterMiniCollection - outerChildAlphasCost ) / totalCost;
    positionGain = ( positionCostOfOuterMiniCollection - outerChildPositionCost ) / totalCost ;


#endif
#endif



    // Return
    if ( std::isnan( miniDataCost + miniAlphasCost + miniPositionCost ) )
    {
      return 0;
    }
    else
    {
      return outerMiniCollection;
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
    std::cout << "\n\n\n\n\n==================================================" << std::endl;
    std::cout << "Trying to expand retained edge " << edgeId << std::endl;

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
  ::OptimizeReferencePosition( AtlasMeshCollection* meshCollection,
                               AtlasMesh::CellIdentifier  vertexId, bool optimize ) const
  {

    // Get the mini collection surrounding the vertex to optimize
#if 0
    AtlasMeshCollection::Pointer  outerMiniCollection = 0;
    AtlasMeshCollection::Pointer  innerMiniCollection = 0;
    if ( !m_AllowMiniMeshSpeedUp )
    {
      // The mini-mesh is the entire mesh
      outerMiniCollection = meshCollection;
      innerMiniCollection = meshCollection;
    }
    else
    {
      // The mini-mesh is an area surrounding the edge to collapse
      innerMiniCollection = meshCollection->GetRegionGrown( vertexId, m_MiniMeshRadius );
      outerMiniCollection = meshCollection->GetRegionGrown( vertexId, m_MiniMeshRadius+1 );
    }
#else
    AtlasMeshCollection::Pointer  outerMiniCollection = meshCollection;
    AtlasMeshCollection::Pointer  innerMiniCollection = meshCollection;
#endif


    // Retrieve the point id of the point corresponding to the vertex
    const AtlasMesh::PointIdentifier  pointId =
      *( outerMiniCollection->GetCells()->ElementAt( vertexId )->PointIdsBegin() );
#if 0
    if ( pointId == 37 )
    {
      outerMiniCollection->Write( "miniBefore.txt" );
      std::cout << "Wrote out miniBefore.txt !!!!!!!!!!!!!!!!!!!!" << std::endl;
    }
#endif

    // Set up the cost calculator
    kvl::AtlasMeshCollectionReferencePositionCost::ParametersType  initialPosition( 3 );
    initialPosition[ 0 ] = outerMiniCollection->GetReferencePosition()->ElementAt( pointId )[ 0 ];
    initialPosition[ 1 ] = outerMiniCollection->GetReferencePosition()->ElementAt( pointId )[ 1 ];
    initialPosition[ 2 ] = outerMiniCollection->GetReferencePosition()->ElementAt( pointId )[ 2 ];
    kvl::AtlasMeshCollectionReferencePositionCost::Pointer  costFunction =
      kvl::AtlasMeshCollectionReferencePositionCost::New();
    costFunction->SetLabelImages( m_LabelImages );
    costFunction->SetPointId( pointId );
    costFunction->SetInitialMeshCollections( innerMiniCollection, outerMiniCollection );
    kvl::AtlasMeshCollectionReferencePositionCost::ParametersType  optimalPosition = initialPosition;

#if 1 /* !!!!!!!!!!!!!!!!! */
    // Optimize the position of the vertex in the reference position
    const bool canMoveX = outerMiniCollection->GetPointParameters()->ElementAt( pointId ).m_CanMoveX;
    const bool canMoveY = outerMiniCollection->GetPointParameters()->ElementAt( pointId ).m_CanMoveY;
    const bool canMoveZ = outerMiniCollection->GetPointParameters()->ElementAt( pointId ).m_CanMoveZ;
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

      std::cout << "parameterOrder: " << parameterOrder << std::endl;


      // Optimize the reference position
      std::cout << "Optimizing the position of the unified point " << pointId << std::endl;
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
        std::cout << "                               Changed reference position for the point " << pointId
                  << " from " << initialPosition << " to " << optimalPosition << std::endl;
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
    //if ( !costFunction->GetValue( optimalPosition ) )
    if ( costFunction->GetValue( optimalPosition ) == itk::NumericTraits< float >::max() )
    {
      return false;
    }
    AtlasMeshCollection::ConstPointer  optimizedOuterCollection = costFunction->GetCostCalculationMeshCollection();
#if 0
    if ( pointId == 37 )
    {
      optimizedOuterCollection->Write( "miniAfter.txt" );
      std::cout << "Wrote out mini after !!!!!!!!!!!!!!!!!!!!" << std::endl;
    }
#endif
    //{
//   float dataCost;
//   float alphasCost;
//   this->GetDataCostAndAlphasCost( optimizedOuterCollection, dataCost, alphasCost );
//   float positionCost = this->GetPositionCost( optimizedOuterCollection );
//   std::cout << "                        Retreiving optimized mesh. Cost is now "
//             << dataCost + alphasCost + positionCost << std::endl;
//   }

    // Copy the positions estimated in the mini mesh to the whole collection
    std::vector< AtlasMesh::PointsContainer::Pointer >  miniPositions =  optimizedOuterCollection->GetPositions();
    std::cout << "Copying positions from mini mesh collection" << std::endl;
    for ( unsigned int meshNumber = 0; meshNumber < miniPositions.size(); meshNumber++ )
    {
      for ( AtlasMesh::PointsContainer::ConstIterator  positionIt = miniPositions[ meshNumber ]->Begin();
            positionIt !=  miniPositions[ meshNumber ]->End(); ++positionIt )
      {
        meshCollection->GetPositions()[ meshNumber ]->ElementAt( positionIt.Index() ) =  positionIt.Value();
      }

    }



    // Also copy the reference position from the mini mesh to the original mesh
    //std::cout << "Copying the reference position from the mini-mesh to the result mesh" << std::endl;
    meshCollection->GetReferencePosition()->ElementAt( pointId ) =
      optimizedOuterCollection->GetReferencePosition()->ElementAt( pointId );

    // Force the triangle area parameters to be updated!!!
    meshCollection->SetK( meshCollection->GetK() );


    // Copy the alphas from the mini mesh to the whole mesh
    for ( AtlasMesh::PointDataContainer::ConstIterator  paramIt = optimizedOuterCollection->GetPointParameters()->Begin();
          paramIt != optimizedOuterCollection->GetPointParameters()->End();
          ++paramIt )
    {
      meshCollection->GetPointParameters()->ElementAt( paramIt.Index() ).m_Alphas = paramIt.Value().m_Alphas;
    }

    return true;
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

      std::cout << "parameterOrder: " << parameterOrder << std::endl;


      // Optimize the reference position
      std::cout << "Optimizing the position of the unified point " << pointId << std::endl;
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
        std::cout << "                               Changed reference position for the point " << pointId
                  << " from " << initialPosition << " to " << optimalPosition << std::endl;
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



#if 0
//
//
//
  AtlasMeshCollection::ConstPointer
  AtlasMeshBuilder
  ::OptimizeForK( const AtlasMeshCollection*  initialMeshCollection ) const
  {
    // Set up the cost calculator
    AtlasMeshCollectionKCost::ParametersType  initialLogK( 1 );
    initialLogK[ 0 ] = log( initialMeshCollection->GetK() );
    AtlasMeshCollectionKCost::Pointer  costFunction = AtlasMeshCollectionKCost::New();
    costFunction->SetLabelImages( m_LabelImages );
    costFunction->SetInitialMeshCollection( initialMeshCollection );
    AtlasMeshCollectionKCost::ParametersType  optimalLogK = initialLogK;

    // Optimize for K
    std::cout << "Optimizing for K" << std::endl;
    itk::PowellOptimizer::Pointer  optimizer = itk::PowellOptimizer::New();
    optimizer->SetCostFunction( costFunction );
    optimizer->SetInitialPosition( initialLogK );
    optimizer->SetAbsolutePrecisionBrent( 0.1 );
    optimizer->SetAbsolutePrecision( 0.1 );
    optimizer->AcceptNewDirectionsOff();
    optimizer->StartOptimization();

    // Retrieve the optimal reference position
    if ( optimizer->GetCurrentValue() != itk::NumericTraits< float >::max() )
    {
      optimalLogK = optimizer->GetCurrentPosition();
      std::cout << "Changed K from " << exp( initialLogK[ 0 ] ) << " to " << exp( optimalLogK[ 0 ] ) << std::endl;
      return costFunction->GetOptimizedMeshCollection( optimalLogK );
    }
    else
    {
      std::cout << "Optimization of K dramatically failed!" << std::endl;
      return initialMeshCollection;
    }


  }
#endif



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
  void
  AtlasMeshBuilder
  ::AnalyzeEdge( AtlasMesh::CellIdentifier edgeId )
  {

#if 1
    mutex.Lock();
#endif

    // Check if the edge is still present in the mesh.
    if ( !m_Current->GetCells()->IndexExists( edgeId ) )
    {
      return;
    }


    std::cout << "    Analyzing edge with id: " << edgeId << std::endl;


    // Get the mini collection surrounding the edge to try to collapse
    AtlasMeshCollection::ConstPointer  outerMiniCollection = m_Current->GetRegionGrown( edgeId, 2 ).GetPointer();
    AtlasMeshCollection::ConstPointer  innerMiniCollection = m_Current->GetRegionGrown( edgeId, 1 ).GetPointer();

#if 1
    mutex.Unlock();
#endif

    //  Calculate cost when you just optimize this edge's vertex positions
    float  retainedCost = 0;
    float  retainedDataCost = 0;
    float  retainedAlphasCost = 0;
    float  retainedPositionCost = 0;
    AtlasMeshCollection::Pointer  retainedMiniCollection =
      this->TryToRetain( innerMiniCollection, outerMiniCollection, edgeId,
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
    AtlasMeshCollection::Pointer  collapsedMiniCollection =
      this->TryToCollapse( innerMiniCollection, outerMiniCollection, edgeId,
                           collapsedDataCost, collapsedAlphasCost, collapsedPositionCost );
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
    std::cout << "                 retainedCost : " << retainedCost
              << "  (" << retainedDataCost << " + " << retainedAlphasCost
              << " + " << retainedPositionCost <<") " << std::endl;
    std::cout << "                 collapsedCost: " << collapsedCost
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

#if 1
    if ( minTotalCostIndex == -1 )
    {
#if 1
      mutex.Lock();
#endif

      std::cout << "Impossible configuration encountered at eget with id " << edgeId << std::endl;
      std::ostringstream  impossibleStream;
      impossibleStream << "impossible_" << edgeId;
      m_Current->Write( impossibleStream.str().c_str() );
      //minTotalCostIndex = 0;

#if 1
      mutex.Unlock();
#endif

      return;
    }

#endif


#if 1
    mutex.Lock();
#endif

    // Do the best move
    if ( minTotalCostIndex == 0 )
    {
      std::cout << "        => retaining edge is best solution" << std::endl;

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
      std::cout << "        => collapsing edge is best solution" << std::endl;

      // std::set< AtlasMesh::CellIdentifier >  disappearingCells;
      // AtlasMesh::CellIdentifier  unifiedVertexIdDummy;
      // this->ExpandCollapse( m_Current, m_EdgeId, m_CollapsedMiniCollection,
      //                         newCurrent, &disappearingCells, unifiedVertexIdDummy );

      // Look up the ids of the two points on the edge
      AtlasMesh::CellType::PointIdConstIterator  pointIt = m_Current->GetCells()->ElementAt( edgeId )->PointIdsBegin();
      const AtlasMesh::PointIdentifier  edgePoint0Id = *pointIt;
      ++pointIt;
      const AtlasMesh::PointIdentifier  edgePoint1Id = *pointIt;


      // What follows will destroy the m_CollapsedMiniCollection, but we need to make sure we calculate its
      // cell links before that happpens
      AtlasMesh::CellLinksContainerPointer  collapsedMiniCollectionCellLinks = collapsedMiniCollection->GetCellLinks();


      // Do the basic mesh surgery: loop over all cells in the original outer mini mesh, and look for them
      // in the collapsed outer mini mesh. If it is found there, delete the corresponding cell in the global mesh,
      // put its pointer to the cell in the collapsed mini mesh, and remove the cell from the collapsed mini mesh
      // (the latter step is necessary as the destructor of the collapsed mini mesh will call delete on the pointer).
      // If it is not found, the cell has disappeared due to the edge collapse, so delete the corresponding cell from
      // the global mesh
      for ( AtlasMesh::CellsContainer::ConstIterator  cellIt = outerMiniCollection->GetCells()->Begin();
            cellIt != outerMiniCollection->GetCells()->End(); ++cellIt )
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
      for ( AtlasMesh::PointsContainer::ConstIterator pointIt = innerMiniCollection->GetReferencePosition()->Begin();
            pointIt != innerMiniCollection->GetReferencePosition()->End(); ++ pointIt )
      {
        m_Current->GetCellLinks()->ElementAt( pointIt.Index() ) =
          collapsedMiniCollectionCellLinks->ElementAt( pointIt.Index() );

        // std::cout << pointIt.Index() << ": copied "
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

#if 1
    mutex.Unlock();
#endif

  }



//
//
//
  bool
  AtlasMeshBuilder
  ::LoadBalancedAnalyzeEdge( std::set< AtlasMesh::CellIdentifier >&  edges,
                             std::map< AtlasMesh::PointIdentifier, int >&  pointOccupancies,
                             int  threadId )
  {

#if 1
    mutex.Lock();
#endif

    // Check if there is anything left to do
    if ( edges.size() == 0 )
    {
#if 1
      mutex.Unlock();
#endif
      return false;
    }


    // Look for the first edge that can be safely worked on
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
#if 1
      mutex.Unlock();
#endif
      sleep( 1 );
      return true;
    }


    AtlasMesh::CellType::PointIdConstIterator  pit
    = m_Current->GetCells()->ElementAt( edgeId )->PointIdsBegin();
    const AtlasMesh::PointIdentifier  p0Id = *pit;
    ++pit;
    const AtlasMesh::PointIdentifier  p1Id = *pit;
    std::cout << "    [THREAD " << threadId << "] " << "    Analyzing edge with id: " << edgeId
              << " (pointIds " << p0Id << " and " << p1Id << " ) " << std::endl;


    // Remove this edge from the edges to be analyzed
    std::cout << "    [THREAD " << threadId << "] " << "         Removing edge with id: " << edgeId << " from edges (" << edges.size() << ")" << std::endl;
    edges.erase( edgeId );
    std::cout << "    [THREAD " << threadId << "] " << "         Removed edge with id: " << edgeId << " from edges (" << edges.size() << ")" << std::endl;


    // Get the mini collection surrounding the edge to try to collapse
    AtlasMeshCollection::ConstPointer  outerMiniCollection = m_Current->GetRegionGrown( edgeId, 2 ).GetPointer();
    AtlasMeshCollection::ConstPointer  innerMiniCollection = m_Current->GetRegionGrown( edgeId, 1 ).GetPointer();


    // Collect the points that are affected by our edge
    std::vector< AtlasMesh::PointIdentifier >  affectedPoints;
    for ( AtlasMesh::PointsContainer::ConstIterator  refPosIt = outerMiniCollection->GetReferencePosition()->Begin();
          refPosIt != outerMiniCollection->GetReferencePosition()->End(); ++refPosIt )
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
#endif


#if 0
    AtlasMeshCollectionValidator::Pointer  validator = AtlasMeshCollectionValidator::New();
    if ( !validator->Validate( outerMiniCollection ) )
    {
      std::ostringstream  outerMiniCollectionStream;
      outerMiniCollectionStream << "outerMiniCollectionStream_thread" << threadId << "_edgeId"<< edgeId << ".txt";
      std::cout << "    [THREAD " << threadId << "] " << "Writing out " << outerMiniCollectionStream.str() << std::endl;
      outerMiniCollection->Write( outerMiniCollectionStream.str().c_str() );
      std::cout << "    [THREAD " << threadId << "] " << "Wrote out " << outerMiniCollectionStream.str() << std::endl;
      exit( -1 );
    }
#endif


    //  Calculate cost when you just optimize this edge's vertex positions
    float  retainedCost = 0;
    float  retainedDataCost = 0;
    float  retainedAlphasCost = 0;
    float  retainedPositionCost = 0;
    AtlasMeshCollection::Pointer  retainedMiniCollection =
      this->TryToRetain( innerMiniCollection, outerMiniCollection, edgeId,
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
    AtlasMeshCollection::Pointer  collapsedMiniCollection =
      this->TryToCollapse( innerMiniCollection, outerMiniCollection, edgeId,
                           collapsedDataCost, collapsedAlphasCost, collapsedPositionCost );
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

#if 1
    if ( minTotalCostIndex == -1 )
    {
#if 1
      mutex.Lock();
#endif

      std::cout << "    [THREAD " << threadId << "] " << "Impossible configuration encountered at eget with id " << edgeId << std::endl;
      std::ostringstream  impossibleStream;
      impossibleStream << "impossible_" << edgeId;
      m_Current->Write( impossibleStream.str().c_str() );
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
#endif

      return true;
    }

#endif


#if 1
    mutex.Lock();
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


      // What follows will destroy the m_CollapsedMiniCollection, but we need to make sure we calculate its
      // cell links before that happpens
      AtlasMesh::CellLinksContainerPointer  collapsedMiniCollectionCellLinks = collapsedMiniCollection->GetCellLinks();


      // Do the basic mesh surgery: loop over all cells in the original outer mini mesh, and look for them
      // in the collapsed outer mini mesh. If it is found there, delete the corresponding cell in the global mesh,
      // put its pointer to the cell in the collapsed mini mesh, and remove the cell from the collapsed mini mesh
      // (the latter step is necessary as the destructor of the collapsed mini mesh will call delete on the pointer).
      // If it is not found, the cell has disappeared due to the edge collapse, so delete the corresponding cell from
      // the global mesh
      for ( AtlasMesh::CellsContainer::ConstIterator  cellIt = outerMiniCollection->GetCells()->Begin();
            cellIt != outerMiniCollection->GetCells()->End(); ++cellIt )
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
      for ( AtlasMesh::PointsContainer::ConstIterator pointIt = innerMiniCollection->GetReferencePosition()->Begin();
            pointIt != innerMiniCollection->GetReferencePosition()->End(); ++ pointIt )
      {
        m_Current->GetCellLinks()->ElementAt( pointIt.Index() ) =
          collapsedMiniCollectionCellLinks->ElementAt( pointIt.Index() );

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
    if ( !validator->Validate( m_Current ) )
    {
      std::cout << "    [THREAD " << threadId << "] " << " invalidated m_Current! " << std::endl;

      std::ostringstream  outerMiniCollectionStream;
      outerMiniCollectionStream << "outerMiniCollectionStream_thread" << threadId << "_edgeId"<< edgeId << ".txt";
      std::cout << "    [THREAD " << threadId << "] " << "Writing out " << outerMiniCollectionStream.str() << std::endl;
      outerMiniCollection->Write( outerMiniCollectionStream.str().c_str() );
      std::cout << "    [THREAD " << threadId << "] " << "Wrote out " << outerMiniCollectionStream.str() << std::endl;

      std::ostringstream  currentStream;
      currentStream << "m_Current_thread" << threadId << "_edgeId"<< edgeId << ".txt";
      std::cout << "    [THREAD " << threadId << "] " << "Writing out " << currentStream.str() << std::endl;
      m_Current->Write( currentStream.str().c_str() );
      std::cout << "    [THREAD " << threadId << "] " << "Wrote out " << currentStream.str() << std::endl;

      exit( -1 );
    }
#endif


    // Now "un-protect" the points we flagged as being worked on
    for ( std::vector< AtlasMesh::PointIdentifier >::const_iterator  it = affectedPoints.begin();
          it != affectedPoints.end(); ++it )
    {
      //std::cout << "    [THREAD " << threadId << "] " << "        Edge " << edgeId << " is now unlocking point " << *it << std::endl;
      ( pointOccupancies[ *it ] )--;
    }


#if 1
    mutex.Unlock();
#endif

    return true;

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
    //mutex.Lock();
    {
      std::ostringstream  descriptionStream;
      descriptionStream << "    [THREAD " << threadId << "] locked mutex because trying to find edge to analyze";
      mutex.DescriptiveLock( descriptionStream.str() );
    }
#endif

    // Check if there is anything left to do
    if ( edges.size() == 0 )
    {
#if 1
      //mutex.Unlock();
      mutex.DescriptiveUnlock();
#endif
      return false;
    }


    // Look for the first edge that can be safely worked on
    std::cout << "    [THREAD " << threadId << "] " << " looking for the first edge that can be safely worked on" << std::endl;
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
#if 1
      //mutex.Unlock();
      mutex.DescriptiveUnlock();
#endif
      sleep( 1 );
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
    //mutex.Unlock();
    mutex.DescriptiveUnlock();
#endif


#if 0
    AtlasMeshCollectionValidator::Pointer  validator = AtlasMeshCollectionValidator::New();
    if ( !validator->Validate( miniCollection ) )
    {
      std::ostringstream  miniCollectionStream;
      miniCollectionStream << "miniCollectionStream_thread" << threadId << "_edgeId"<< edgeId << ".txt";
      std::cout << "    [THREAD " << threadId << "] " << "Writing out " << miniCollectionStream.str() << std::endl;
      miniCollection->Write( miniCollectionStream.str().c_str() );
      std::cout << "    [THREAD " << threadId << "] " << "Wrote out " << miniCollectionStream.str() << std::endl;

      std::ostringstream  currentStream;
      currentStream << "m_Current_thread" << threadId << "_edgeId"<< edgeId << ".txt";
      std::cout << "    [THREAD " << threadId << "] " << "Writing out " << currentStream.str() << std::endl;
      m_Current->Write( currentStream.str().c_str() );
      std::cout << "    [THREAD " << threadId << "] " << "Wrote out " << currentStream.str() << std::endl;

      exit( -1 );
    }
#endif



    itk::TimeProbe  timeProbe;
    timeProbe.Start();


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


    timeProbe.Stop();
    std::cout << "    [THREAD " << threadId << "] took " << timeProbe.GetMeanTime() << " seconds to evaluate moves" << std::endl;


#if 1
    if ( minTotalCostIndex == -1 )
    {
#if 1
      //mutex.Lock();
      {
        std::ostringstream  descriptionStream;
        descriptionStream << "    [THREAD " << threadId << "] locked mutex because impossible configuration encountered";
        mutex.DescriptiveLock( descriptionStream.str() );
      }
#endif

      std::cout << "    [THREAD " << threadId << "] " << "Impossible configuration encountered at eget with id " << edgeId << std::endl;
      std::ostringstream  impossibleStream;
      impossibleStream << "impossible_" << edgeId;
      m_Current->Write( impossibleStream.str().c_str() );
      //minTotalCostIndex = 0;

      // Now "un-protect" the points we flagged as being worked on
      for ( std::vector< AtlasMesh::PointIdentifier >::const_iterator  it = affectedPoints.begin();
            it != affectedPoints.end(); ++it )
      {
        //std::cout << "    [THREAD " << threadId << "] " << "        Edge " << edgeId << " is now unlocking point " << *it << std::endl;
        ( pointOccupancies[ *it ] )--;
      }

#if 1
      //mutex.Unlock();
      mutex.DescriptiveUnlock();
#endif

      return true;
    }

#endif


#if 1
    //mutex.Lock();
    {
      std::ostringstream  descriptionStream;
      descriptionStream << "    [THREAD " << threadId << "] locked mutex because applying best move to the global mesh";
      mutex.DescriptiveLock( descriptionStream.str() );
    }
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

#if 0
    if ( !validator->Validate( m_Current ) )
    {
      std::cout << "    [THREAD " << threadId << "] " << " invalidated m_Current! " << std::endl;

      std::ostringstream  miniCollectionStream;
      miniCollectionStream << "miniCollectionStream_thread" << threadId << "_edgeId"<< edgeId << ".txt";
      std::cout << "    [THREAD " << threadId << "] " << "Writing out " << miniCollectionStream.str() << std::endl;
      miniCollection->Write( miniCollectionStream.str().c_str() );
      std::cout << "    [THREAD " << threadId << "] " << "Wrote out " << miniCollectionStream.str() << std::endl;

      std::ostringstream  currentStream;
      currentStream << "m_Current_thread" << threadId << "_edgeId"<< edgeId << ".txt";
      std::cout << "    [THREAD " << threadId << "] " << "Writing out " << currentStream.str() << std::endl;
      m_Current->Write( currentStream.str().c_str() );
      std::cout << "    [THREAD " << threadId << "] " << "Wrote out " << currentStream.str() << std::endl;

      exit( -1 );
    }
#endif


    // Now "un-protect" the points we flagged as being worked on
    for ( std::vector< AtlasMesh::PointIdentifier >::const_iterator  it = affectedPoints.begin();
          it != affectedPoints.end(); ++it )
    {
      //std::cout << "    [THREAD " << threadId << "] " << "        Edge " << edgeId << " is now unlocking point " << *it << std::endl;
      ( pointOccupancies[ *it ] )--;
    }

    std::cout << "    [THREAD " << threadId << "] " << "Done with edge " << edgeId << std::endl;

#if 1
    //mutex.Unlock();
    mutex.DescriptiveUnlock();
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

    std::cout << "    Analyzing edge with id: " << edgeId << std::endl;


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
    std::cout << "                 retainedCost : " << retainedCost
              << "  (" << retainedDataCost << " + " << retainedAlphasCost
              << " + " << retainedPositionCost <<") " << std::endl;
    std::cout << "                 collapsedCost: " << collapsedCost
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

#if 1
    if ( minTotalCostIndex == -1 )
    {
#if 1
      mutex.Lock();
#endif

      std::cout << "Impossible configuration encountered at eget with id " << edgeId << std::endl;
      std::ostringstream  impossibleStream;
      impossibleStream << "impossible_" << edgeId;
      m_Current->Write( impossibleStream.str().c_str() );
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
      std::cout << "        => retaining edge is best solution" << std::endl;

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
      std::cout << "        => collapsing edge is best solution" << std::endl;

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
    std::cout << "Selecting a subset of edges that can be analyzed independently..." << std::endl;
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

    std::cout << "   resulting number of edges: " << independentEdges.size()
              << " (" << static_cast< float >( independentEdges.size() )
              / static_cast< float >( edgesToTry.size() ) * 100.0f
              << " % of all edges)" << std::endl;


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
  ::ThreaderCallback( void *arg )
  {

    // Retrieve the input arguments
    const int  threadId = ((itk::MultiThreader::ThreadInfoStruct *)(arg))->ThreadID;
    const int  threadCount = ((itk::MultiThreader::ThreadInfoStruct *)(arg))->NumberOfThreads;

    const ThreadStruct*  str = (ThreadStruct *)(((itk::MultiThreader::ThreadInfoStruct *)(arg))->UserData);

#if 0
    //
    std::cout << "I am thread " << threadId << " and I'm analyzing edge " <<  str->m_EdgesToTry[ threadId ] << std::endl;
    str->m_Builder->AnalyzeEdge( str->m_EdgesToTry[ threadId ] );

#else
    //
    for ( unsigned int edgesToTryIndex = threadId; edgesToTryIndex < str->m_EdgesToTry.size(); edgesToTryIndex += threadCount )
    {
      std::cout << "I am thread " << threadId << " and I'm analyzing edge " <<  str->m_EdgesToTry[ edgesToTryIndex ] << std::endl;
      std::cout << "Progress: " << static_cast< float >( edgesToTryIndex )
                / static_cast< float >( str->m_EdgesToTry.size() ) * 100.0f
                << std::endl;

      str->m_Builder->AnalyzeEdge( str->m_EdgesToTry[ edgesToTryIndex ] );
    }
#endif

    return ITK_THREAD_RETURN_VALUE;
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


//
//
//
  ITK_THREAD_RETURN_TYPE
  AtlasMeshBuilder
  ::FastThreaderCallback( void *arg )
  {

    // Retrieve the input arguments
    const int  threadId = ((itk::MultiThreader::ThreadInfoStruct *)(arg))->ThreadID;
    //const int  threadCount = ((itk::MultiThreader::ThreadInfoStruct *)(arg))->NumberOfThreads;

    FastThreadStruct*  str = (FastThreadStruct *)(((itk::MultiThreader::ThreadInfoStruct *)(arg))->UserData);

    // Loop over all edges to analyze, and do it
    std::cout << "I'm thread " << threadId << " and I'm starting to analyze "
              << str->m_Items[ threadId ].m_EdgeIds.size() << " edges" << std::endl;
    for ( unsigned int edgeNumber = 0; edgeNumber < str->m_Items[ threadId ].m_EdgeIds.size(); edgeNumber++ )
    {
      AtlasMeshCollection::ConstPointer  miniCollection = str->m_Items[ threadId ].m_MiniCollections[ edgeNumber ];
      const AtlasMesh::CellIdentifier  edgeId = str->m_Items[ threadId ].m_EdgeIds[ edgeNumber ];
      const AtlasMesh::PointIdentifier  edgePoint0Id = str->m_Items[ threadId ].m_EdgePoint0Ids[ edgeNumber ];
      const AtlasMesh::PointIdentifier  edgePoint1Id = str->m_Items[ threadId ].m_EdgePoint1Ids[ edgeNumber ];

      std::cout << "    thread " << threadId << " analyzing edge: " << edgeId << std::endl;

      str->m_Builder->AnalyzeEdgeFast( miniCollection, edgeId, edgePoint0Id, edgePoint1Id,
                                       str->m_Items[ threadId ].m_DisappearingPointsLookupTable,
                                       str->m_Items[ threadId ].m_NewReferencePosition,
                                       str->m_Items[ threadId ].m_NewPositions,
                                       str->m_Items[ threadId ].m_NewPointParameters,
                                       str->m_Items[ threadId ].m_DisappearingCells );

      std::cout << "Progress: " << static_cast< float >( edgeNumber )
                / static_cast< float >( str->m_Items[ threadId ].m_EdgeIds.size() ) * 100.0f
                << std::endl;
    }

    return ITK_THREAD_RETURN_VALUE;
  }



} // end namespace kvl
