#include "kvlAtlasParameterEstimator.h"

#include "itkImageRegionConstIterator.h"
#include "kvlAtlasMeshLabelImageStatisticsCollector.h"
#include "kvlAtlasMeshToLabelImageCostAndGradientCalculator.h"
#include "kvlAtlasMeshDeformationConjugateGradientOptimizer.h"
#include "kvlAtlasMeshDeformationFixedStepGradientDescentOptimizer.h"
#include "kvlAtlasMeshDeformationGradientDescentOptimizer.h"
#include "kvlAtlasMeshDeformationLBFGSOptimizer.h"
#include "kvlAtlasMeshSmoother.h"
#include "itkCommand.h"


namespace kvl
{


//
//
//
AtlasParameterEstimator
::AtlasParameterEstimator()
{

  m_MeshCollection = 0;
  m_CompressionLookupTable = 0;
  m_IterationNumber = 0;
  m_MaximumNumberOfIterations = 300;
  m_LabelImageNumber = 0;
  m_NumberOfLabelImages = 0;
  m_AlphasEstimationIterationNumber = 0;
  m_AlphasEstimationMaximumNumberOfIterations = 20;
  m_PositionEstimationIterationNumber = 0;
  m_PositionEstimationMaximumNumberOfIterations = 0; // To be retrieved from optimizer
  m_PositionEstimationIterationEventResolution = 10; 
  m_CurrentMinLogLikelihoodTimesPrior = 0.0f;
  m_AlphaEstimationStopCriterion = 0.0005f;
  m_AlphasSmoothingFactor = 0.0f;
  m_StopCriterion = 0.001f;

  m_PositionOptimizer = LBFGS;

  m_NumberOfThreads = itk::MultiThreader::GetGlobalDefaultNumberOfThreads();
  
}



//
//
//
AtlasParameterEstimator
::~AtlasParameterEstimator()
{

}



//
//
//
void 
AtlasParameterEstimator
::PrintSelf( std::ostream& os, itk::Indent indent ) const
{



}



//
//
//
const AtlasParameterEstimator::LabelImageType*
AtlasParameterEstimator
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
AtlasParameterEstimator
::SetLabelImages( const std::vector< LabelImageType::ConstPointer >& labelImages,
                  const CompressionLookupTable*  lookupTable )
{

  m_CompressionLookupTable = lookupTable;
  
  if ( labelImages.size() == 0 )
    return;
  
  m_LabelImages = labelImages;
  m_NumberOfLabelImages = m_LabelImages.size();
 
}




//
//
//
void 
AtlasParameterEstimator
::Estimate( bool verbose )
{
  // Sanity checking
  if ( m_MeshCollection->GetNumberOfMeshes() == 0 )
    {
    itkExceptionMacro( "No initial meshes set!" );
    }    
  if ( m_LabelImages.size() == 0 )
    {
    itkExceptionMacro( "No label images set!" );
    }
  if ( m_MeshCollection->GetNumberOfMeshes() != m_LabelImages.size() )
    {
    itkExceptionMacro( "Number of meshes must match number of label images!" );
    }  
      
  // Main loop  
  double  previousMinLogLikelihoodTimesPrior = 1e15;
  m_CurrentMinLogLikelihoodTimesPrior = previousMinLogLikelihoodTimesPrior / 2;
  m_IterationNumber = 0;
  this->InvokeEvent( itk::StartEvent() );
  while ( ( ( ( previousMinLogLikelihoodTimesPrior - m_CurrentMinLogLikelihoodTimesPrior ) / 
               fabsf( m_CurrentMinLogLikelihoodTimesPrior ) ) > m_StopCriterion ) &&
          ( m_IterationNumber < m_MaximumNumberOfIterations ) ) 
    {
    // Estimate alphas
    this->EstimateAlphas();
 
    // Smooth results
    this->SmoothAlphas();

    // Register. This will also provide us with the minLogLikelihoodTimesPrior of 
    // the current parameter set
    previousMinLogLikelihoodTimesPrior = m_CurrentMinLogLikelihoodTimesPrior;
    m_CurrentMinLogLikelihoodTimesPrior = this->EstimatePositions();
    
    // Prepare for the next iteration
    if ( verbose )
     {
     std::cout << "================ IterationNumber: " << m_IterationNumber 
              << " CurrentMinLogLikelihoodTimesPrior: " << m_CurrentMinLogLikelihoodTimesPrior
              << "================" << std::endl;
      }
    m_IterationNumber++;
    this->InvokeEvent( itk::IterationEvent() );
    }


  // Estimate final alphas
  this->EstimateAlphas();
  this->InvokeEvent( itk::EndEvent() );
  
}



//
//
//
void
AtlasParameterEstimator
::EstimateAlphas()
{

  // If we're smoothing, we have to start with flat alpha entries
  if ( m_AlphasSmoothingFactor != 0 )
    {
    m_MeshCollection->FlattenAlphas();
    }


  // Allocate a container to sum the statistics over all label images
  typedef AtlasMeshLabelImageStatisticsCollector::StatisticsContainerType StatisticsContainerType;
  StatisticsContainerType::Pointer  pooledStatistics = StatisticsContainerType::New();
  AtlasMesh::PointDataContainer::ConstIterator  pointParamIt = m_MeshCollection->GetPointParameters()->Begin();
  while ( pointParamIt != m_MeshCollection->GetPointParameters()->End() )
    {
    pooledStatistics->CreateIndex( pointParamIt.Index() );
    
    ++pointParamIt;
    }
  
  

  double  previousCost = 1e15;
  double  currentCost = previousCost / 2;
  m_AlphasEstimationIterationNumber = 0;
  this->InvokeEvent( AlphasEstimationStartEvent() );
  while ( ( ( ( previousCost - currentCost ) / fabsf( currentCost ) ) > m_AlphaEstimationStopCriterion ) && 
          ( m_AlphasEstimationIterationNumber < m_AlphasEstimationMaximumNumberOfIterations ) ) 
    {
    previousCost = currentCost;
    
    // Initialize pooled statistics to zero
    unsigned int  numberOfLabeles = m_MeshCollection->GetPointParameters()->Begin().Value().m_Alphas.Size();
    StatisticsContainerType::Element  zeroEntry( numberOfLabeles );
    zeroEntry.Fill( 0.0f );
    StatisticsContainerType::Iterator  pooledIt = pooledStatistics->Begin();
    while ( pooledIt != pooledStatistics->End() )
      {
      pooledIt.Value() = zeroEntry;
      pooledIt++;
      }


    // Loop over all label images, retrieve their statistics, and add those to the pooled statistics
    currentCost = 0;
    int  numberOfLabelImagesToVisit = m_NumberOfLabelImages;
    AtlasMeshLabelImageStatisticsCollector::Pointer  statisticsCollector = AtlasMeshLabelImageStatisticsCollector::New(); 
    statisticsCollector->SetNumberOfThreads( m_NumberOfThreads );
    for ( int  labelImageNumber = 0; labelImageNumber < numberOfLabelImagesToVisit; labelImageNumber++ )
      {
      // Calculate statistics for this label image 
      statisticsCollector->SetLabelImage( m_LabelImages[ labelImageNumber ], m_CompressionLookupTable );
      statisticsCollector->Rasterize( m_MeshCollection->GetMesh( labelImageNumber ) );

      // Add statistics to pooled statistics
      pooledIt = pooledStatistics->Begin();
      StatisticsContainerType::ConstIterator  statIt = statisticsCollector->GetLabelStatistics()->Begin();
      while ( pooledIt != pooledStatistics->End() )
        {
        pooledIt.Value() += statIt.Value();

        pooledIt++;
        statIt++;
        }
        
      // Also add contribution of cost function
      currentCost += statisticsCollector->GetMinLogLikelihood();
      }


    // Estimated alphas are simply normalized pooled statistics
    pooledIt = pooledStatistics->Begin();
    PointDataContainerType::Iterator  pointParamIt = m_MeshCollection->GetPointParameters()->Begin();
    while ( pooledIt != pooledStatistics->End() )
      {
      if ( pointParamIt.Value().m_CanChangeAlphas )
        {
        pointParamIt.Value().m_Alphas = pooledIt.Value();
        pointParamIt.Value().m_Alphas /= pooledIt.Value().sum();
        }
        
      pooledIt++;
      pointParamIt++;
      }

     
    // Prepare for next iteration
    m_AlphasEstimationIterationNumber++;
    this->InvokeEvent( AlphasEstimationIterationEvent() );
    
    } // End iterative loop
  
  this->InvokeEvent( AlphasEstimationEndEvent() );
    

}
    
  

  
  
//
//
//
void
AtlasParameterEstimator
::SmoothAlphas()
{

  if ( m_AlphasSmoothingFactor != 0 )
    {
    // Retrieve sigma of the smoothing kernel
    const double  maximalSigma = 0.05 * ( m_LabelImages[ 0 ]->GetLargestPossibleRegion().GetSize()[ 0 ] );
    const double  sigma = maximalSigma * ( 1 - log( 2 - m_AlphasSmoothingFactor ) / log( 2 ) );
    
    // Smooth
    AtlasMeshSmoother::Pointer  smoother = AtlasMeshSmoother::New();
    smoother->SetMeshCollection( m_MeshCollection );
    smoother->SetSigma( sigma );
    m_MeshCollection = smoother->GetSmoothedMeshCollection();

#if 1
    {
    // std::cout << "sigma: " << sigma << std::endl;
    std::ostringstream  fileNameStream;
    fileNameStream << "smoothed_" << sigma << ".txt";
    m_MeshCollection->Write( fileNameStream.str().c_str() ); 
    }
#endif
    }

}




//
//
//  
double
AtlasParameterEstimator
::EstimatePositions()
{

  double  totalMinLogLikelihoodTimesPrior = 0.0f;
  
  for ( m_LabelImageNumber = 0; m_LabelImageNumber < m_NumberOfLabelImages; m_LabelImageNumber++ )
    {
    //
    totalMinLogLikelihoodTimesPrior += this->EstimatePosition( m_LabelImageNumber );
    }

  return totalMinLogLikelihoodTimesPrior;
}



//
//
//  
double
AtlasParameterEstimator
::EstimatePosition( unsigned int  labelImageNumber )
{

  // Let's check if we are expected to just stay put
  if ( m_MeshCollection->GetK() >= 1000 )
    {
    // Just return the cost at this position
    double  cost;
    this->CalculateCurrentPositionCostAndGradient( labelImageNumber, cost );
    return cost;
    }


  // Set up gradient calculator
  AtlasMeshToLabelImageCostAndGradientCalculator::Pointer  
              calculator = AtlasMeshToLabelImageCostAndGradientCalculator::New();  
  calculator->SetLabelImage( m_LabelImages[ labelImageNumber ] , 
                             m_CompressionLookupTable );
  calculator->SetNumberOfThreads( m_NumberOfThreads );
  
  
  // Set up a deformation optimizer and pass the gradient calculator on to it
  AtlasMeshDeformationOptimizer::Pointer  optimizer = 0;
  switch( m_PositionOptimizer ) 
    {
    case FIXED_STEP_GRADIENT_DESCENT: 
      {
      optimizer = AtlasMeshDeformationFixedStepGradientDescentOptimizer::New();
      optimizer->SetVerbose( false );
      break;
      } 
    case GRADIENT_DESCENT: 
      {
      optimizer = AtlasMeshDeformationGradientDescentOptimizer::New();
      optimizer->SetVerbose( true );
      break;
      } 
    case CONJUGATE_GRADIENT: 
      {
      optimizer = AtlasMeshDeformationConjugateGradientOptimizer::New();
      optimizer->SetVerbose( true );
      break;
      } 
    default:
      {
      optimizer = AtlasMeshDeformationLBFGSOptimizer::New();
      optimizer->SetVerbose( true );
      break;
      }
    }
  optimizer->SetCostAndGradientCalculator( calculator );
  optimizer->SetMesh( const_cast< AtlasMesh* >( m_MeshCollection->GetMesh( labelImageNumber ) ) );
  
  // Also make sure iteration events generated by the optimizer are forwarded to our own users
  typedef itk::MemberCommand< Self >   MemberCommandType;
  MemberCommandType::Pointer  command = MemberCommandType::New();
  command->SetCallbackFunction( this, &AtlasParameterEstimator::HandleOptimizerEvent );
  optimizer->AddObserver( DeformationStartEvent(), command );
  optimizer->AddObserver( DeformationIterationEvent(), command );
  optimizer->AddObserver( DeformationEndEvent(), command );
  m_PositionEstimationMaximumNumberOfIterations = optimizer->GetMaximumNumberOfIterations();
  optimizer->SetIterationEventResolution( m_PositionEstimationIterationEventResolution );
  
  // Start the optimization
  optimizer->Go();
  
  // We've given the mesh a new position container, but the meshCollection doesn't know about this!
  m_MeshCollection->GetPositions()[ labelImageNumber ] = const_cast< AtlasMesh::PointsContainer* >( optimizer->GetMesh()->GetPoints() );
    
  return optimizer->GetMinLogLikelihoodTimesPrior();

}



  

//
//
//
AtlasPositionGradientContainerType::Pointer  
AtlasParameterEstimator
::CalculateCurrentPositionCostAndGradient( unsigned int labelImageNumber, double& minLogLikelihoodTimesPrior ) const
{
  
  // Sanity checking
  if ( ( labelImageNumber >= m_MeshCollection->GetNumberOfMeshes() ) || 
       ( labelImageNumber >= m_LabelImages.size() ) )
    {
    return 0;
    }    

  // Set up gradient calculator
  AtlasMeshToLabelImageCostAndGradientCalculator::Pointer  
              calculator = AtlasMeshToLabelImageCostAndGradientCalculator::New();  
  calculator->SetLabelImage( m_LabelImages[ labelImageNumber ] , 
                             m_CompressionLookupTable );
  calculator->SetNumberOfThreads( m_NumberOfThreads );
  
  // Let it do its work
  calculator->Rasterize( m_MeshCollection->GetMesh( labelImageNumber ) );
  minLogLikelihoodTimesPrior = calculator->GetMinLogLikelihoodTimesPrior();

  return calculator->GetPositionGradient();

}



//
//
//
AtlasPositionGradientContainerType::Pointer
AtlasParameterEstimator
::GetCurrentPositionGradient( unsigned int labelImageNumber ) const
{
  double  dummy;
  
  return this->CalculateCurrentPositionCostAndGradient( labelImageNumber, dummy );
}



//
//
//
void
AtlasParameterEstimator
::HandleOptimizerEvent( itk::Object* object, const itk::EventObject & event )
{
  
  if ( typeid( event ) == typeid( DeformationStartEvent ) )
    {
    m_PositionEstimationIterationNumber = 0;
    this->InvokeEvent( PositionEstimationStartEvent() );
    }
  else if ( typeid( event ) == typeid( DeformationIterationEvent ) )
    {
    this->InvokeEvent( PositionEstimationIterationEvent() );
    m_PositionEstimationIterationNumber +=  m_PositionEstimationIterationEventResolution;
    }
  else if ( typeid( event ) == typeid( DeformationEndEvent ) )
    {
    this->InvokeEvent( PositionEstimationEndEvent() );
    }

  
}   
   

} // end namespace kvl

