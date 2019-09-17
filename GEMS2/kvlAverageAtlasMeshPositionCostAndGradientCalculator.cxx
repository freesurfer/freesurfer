#include "kvlAverageAtlasMeshPositionCostAndGradientCalculator.h"

#include "kvlAtlasMeshCollection.h"


namespace kvl
{

//
//
//
AverageAtlasMeshPositionCostAndGradientCalculator
::AverageAtlasMeshPositionCostAndGradientCalculator()
{
}


//
//
//
AverageAtlasMeshPositionCostAndGradientCalculator
::~AverageAtlasMeshPositionCostAndGradientCalculator()
{
}




//
//
//
void 
AverageAtlasMeshPositionCostAndGradientCalculator
::Rasterize( const AtlasMesh* mesh )
{

  //
  //std::cout << "Starting rasterization" << std::endl;
  
  // Sanity check
  if ( ( m_Positions.size() == 0 ) ||
       ( m_Ks.size() != m_Positions.size() ) )
    {
    itkExceptionMacro( << "No positions set or Ks don't not same number as positions" );  
    }  
  
  // Initialize from a clean slate
  //std::cout << "Initializing" << std::endl;
  m_MinLogLikelihoodTimesPrior = 0;
  m_PositionGradient = AtlasPositionGradientContainerType::New();
  AtlasPositionGradientType  zeroEntry( 0.0f );
  for ( AtlasMesh::PointsContainer::ConstIterator pointIt = mesh->GetPoints()->Begin();
        pointIt != mesh->GetPoints()->End(); ++pointIt )
    {
    m_PositionGradient->InsertElement( pointIt.Index(), zeroEntry );
    }


  // For each positon/K pair, generate a mesh of which we're only ever going to change 
  // a pointer to its position container. For efficiency purposes we'll cache these since
  // computing each mesh's ReferenceTetrahedronInfo's (cellData) takes some time
  if ( m_CachedInternalMeshes.size() == 0 )
    {
    //std::cout << "Building internal meshes" << std::endl;
    for ( int meshNumber = 0; meshNumber < m_Positions.size(); meshNumber++ )
      {
      //std::cout << "  meshNumber: " << meshNumber << std::endl;  
      AtlasMeshCollection::Pointer  helper = AtlasMeshCollection::New(); 
      helper->GenerateFromSingleMesh( const_cast< AtlasMesh* >( mesh ), 0, m_Ks[ meshNumber ] );
      helper->SetReferencePosition( const_cast< AtlasMesh::PointsContainer* >( m_Positions[ meshNumber ].GetPointer() ) );
      m_CachedInternalMeshes.push_back( const_cast< AtlasMesh* >( helper->GetReferenceMesh().GetPointer() ) );
      }
    }    

  // Now loop over all meshes, each time setting the mesh's position to the current one, and then computing the cost
  // and gradient of only the deformation prior
  //std::cout << "Rasterizing internal meshes" << std::endl;
  for ( int meshNumber = 0; meshNumber < m_CachedInternalMeshes.size(); meshNumber++ )
    {
    //std::cout << "  meshNumber: " << meshNumber << std::endl; 
  
    // Set up mesh to rasterize
    AtlasMesh::Pointer  internalMesh = m_CachedInternalMeshes[ meshNumber ];
    internalMesh->SetPoints( const_cast< AtlasMesh::PointsContainer* >( mesh->GetPoints() ) );
  
    // Now do the work
    AtlasMeshPositionCostAndGradientCalculator::Pointer  calculator = AtlasMeshPositionCostAndGradientCalculator::New();
    calculator->Rasterize( internalMesh );
    
    // If successful, add contribution to cost
    const double  minLogLikelihoodContribution = calculator->GetMinLogLikelihoodTimesPrior();
    if ( minLogLikelihoodContribution == itk::NumericTraits< double >::max() )
      {
      // Problem with this mesh - abort
      //std::cout << "Aborting" << std::endl;
      m_MinLogLikelihoodTimesPrior = itk::NumericTraits< double >::max();
      return;
      }
    else
      {
      m_MinLogLikelihoodTimesPrior += minLogLikelihoodContribution;  
      }
      
    // Also add contribution to gradient  
    AtlasPositionGradientContainerType::ConstIterator  sourceIt = calculator->GetPositionGradient()->Begin(); 
    AtlasPositionGradientContainerType::Iterator  targetIt = m_PositionGradient->Begin(); 
    for ( ; targetIt != m_PositionGradient->End(); ++sourceIt, ++targetIt )
      {
      targetIt.Value() += sourceIt.Value();
      } 
    
    } // End loop over meshes
    
    
  // Make sure boundary conditions are respected
  this->ImposeBoundaryCondition( mesh );
  
  
}



} // end namespace kvl
