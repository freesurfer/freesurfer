#include "kvlAtlasMeshDeformationLevenbergMarquardtOptimizer.h"


namespace kvl
{


//
//
//
AtlasMeshDeformationLevenbergMarquardtOptimizer
::AtlasMeshDeformationLevenbergMarquardtOptimizer()
{
  m_LevenbergMarquardt = AtlasMeshDeformationLevenbergMarquardt::New();
  m_TrialLevenbergMarquardt = AtlasMeshDeformationLevenbergMarquardt::New();

  m_Lambda = 10.0f;
  
  m_MaximalDeformationStopCriterion = 0.05;

}


//
//
//
AtlasMeshDeformationLevenbergMarquardtOptimizer
::~AtlasMeshDeformationLevenbergMarquardtOptimizer()
{
}



//
//
//
void
AtlasMeshDeformationLevenbergMarquardtOptimizer
::Initialize()
{
  
  //std::cout << "Starting initializing" << std::endl;
  m_Lambda = 10.0f;

  // Create an empty, dummy template image needed by anything that rasterizes meshes
  typedef AtlasMeshDeformationLevenbergMarquardt::LabelImageType  LabelImageType;
  LabelImageType::Pointer  dummyTemplateImage = LabelImageType::New();
  if ( m_SegmentedImage )
    {
    dummyTemplateImage->SetRegions( m_SegmentedImage->GetBufferedRegion() );  
    }
  else if ( m_ProbabilityImage )
    {
    dummyTemplateImage->SetRegions( m_ProbabilityImage->GetBufferedRegion() );
    }
  else
    {
    dummyTemplateImage->SetRegions( m_Images[0]->GetBufferedRegion() );
    }
  //dummyTemplateImage->Allocate();

  // Set up m_LevenbergMarquardt
  m_LevenbergMarquardt->SetMapCompToComp(m_mapCompToComp);
  m_LevenbergMarquardt->SetLabelImage( dummyTemplateImage );
  if ( m_SegmentedImage )
    {
    m_LevenbergMarquardt->SetSegmentedImage( m_SegmentedImage );
    }
  else if ( m_ProbabilityImage )
    {
    m_LevenbergMarquardt->SetProbabilityImage( m_ProbabilityImage );
    }
  else
    {
    m_LevenbergMarquardt->SetImages( m_Images );
    m_LevenbergMarquardt->SetMeans( m_Means );
    m_LevenbergMarquardt->SetPrecisions( m_Precisions );
    }
  m_LevenbergMarquardt->SetMeshToImageTransform( m_MeshToImageTransform );
  //m_LevenbergMarquardt->SetMaximumNumberOfCGIterations( 2 );
  
  // Idem for m_TrialLevenbergMarquardt
  m_TrialLevenbergMarquardt->SetLabelImage( dummyTemplateImage );
  if ( m_SegmentedImage )
    {
    m_TrialLevenbergMarquardt->SetSegmentedImage( m_SegmentedImage );  
    }
  else if ( m_ProbabilityImage )
    {
    m_TrialLevenbergMarquardt->SetProbabilityImage( m_ProbabilityImage );
    }
  else
    {
    m_TrialLevenbergMarquardt->SetImages( m_Images );
    m_TrialLevenbergMarquardt->SetMeans( m_Means );
    m_TrialLevenbergMarquardt->SetPrecisions( m_Precisions );
    }
  m_TrialLevenbergMarquardt->SetMeshToImageTransform( m_MeshToImageTransform );
  //m_TrialLevenbergMarquardt->SetMaximumNumberOfCGIterations( 2 );

  // Rasterize m_LevenbergMarquardt
  m_LevenbergMarquardt->Rasterize( m_Mesh );
  
  //m_Lambda = 1.0f;
  m_Initialized = true;
  
  //std::cout << "Done initializing" << std::endl;

}



//
//
//
bool
AtlasMeshDeformationLevenbergMarquardtOptimizer
::Go()
{
  // std::cout << "m_IterationEventResolution: " << m_IterationEventResolution << std::endl;

  this->InvokeEvent( DeformationStartEvent() );

  //std::cout << "m_PositionUpdatingMaximumNumberOfIterations: " << m_PositionUpdatingMaximumNumberOfIterations << std::endl;

  
  bool haveMoved = false; // Keep track if we've ever moved or not
  for ( m_IterationNumber = 0; 
        m_IterationNumber < m_MaximumNumberOfIterations;
        m_IterationNumber++ )
    {
    // Try to make one successful move
    const float maximalDeformation = this->PerformOneSuccessfulStep();
    if ( maximalDeformation > 0 )
      {
      haveMoved = true;
      }
      
    // Test if we need to stop
    if ( maximalDeformation <= m_MaximalDeformationStopCriterion )
      {
//      std::cout << "        maximalDeformation is too small; stopping" << std::endl;
      break;
      }

    if ( !( m_IterationNumber % m_IterationEventResolution ) )
      {
//      std::cout << "Generating DeformationIterationEvent" << std::endl;
      this->InvokeEvent( DeformationIterationEvent() );
      }

    } // End loop over all positive moves


  this->InvokeEvent( DeformationEndEvent() );

  return haveMoved;
 
}



//
//
//
float
AtlasMeshDeformationLevenbergMarquardtOptimizer
::PerformOneSuccessfulStep()
{

#if 0
  m_TrialLevenbergMarquardt->SetPeerToCopyHessianFrom( m_LevenbergMarquardt );
#endif  
  
  const float  lambdaIncreaseFactor = 2.0f; // Factor by which lambda is increased if failure detected
  const float  lambdaDecreaseFactor = 1.1f; // Factor by which lambda is decreased if success detected

  //std::cout << "Starting AtlasMeshDeformationLevenbergMarquardtOptimizer::PerformOneSuccessfulStep" << std::endl;

  
  // Make sure we're initialized
  if ( !m_Initialized )
    {
    this->Initialize();
    }

  // Evaluate current cost
  double  currentCost = m_LevenbergMarquardt->GetMinLogLikelihoodTimesPrior();

  
  // Keep trying until we've found something sucessful (or we give up)
  while ( true )
    {
    // Save current position
    AtlasMesh::PointsContainer::Pointer  currentPosition = m_Mesh->GetPoints();

    // Calculate the trial step with the current lambda
    AtlasPositionGradientContainerType::Pointer  trialStep = m_LevenbergMarquardt->GetStep( m_Lambda, false );

    // Add the step to the current position to obtain the trial position
    float  maximalDeformation = 0.0f;
    AtlasMesh::PointsContainer::Pointer  trialPosition = AtlasMesh::PointsContainer::New();
    AtlasMesh::PointsContainer::ConstIterator  posIt = currentPosition->Begin();
    AtlasPositionGradientContainerType::ConstIterator  stepIt = trialStep->Begin();
    while ( stepIt != trialStep->End() )
      {
      AtlasPositionGradientType  trialStep = stepIt.Value();

      if ( trialStep.GetNorm() > maximalDeformation )
        {
        maximalDeformation = trialStep.GetNorm();
        }

      trialPosition->InsertElement( posIt.Index(), posIt.Value() + trialStep );

      // Move on to next vertex
      ++posIt;
      ++stepIt;
      }

    //std::cout << "maximalDeformation: " << maximalDeformation << std::endl;


    
    
    // Evaluate the trial position
    m_Mesh->SetPoints( trialPosition );
    m_TrialLevenbergMarquardt->SetMapCompToComp(m_mapCompToComp);
    m_TrialLevenbergMarquardt->Rasterize( m_Mesh );
    double  trialCost = m_TrialLevenbergMarquardt->GetMinLogLikelihoodTimesPrior();


    // Do the Right Thing
//    std::cout << "=======================" << std::endl;
    if ( trialCost < currentCost )
      {
      // Good move. Keep the move and decrease Marquardt's lambda
//      std::cout << "Good trial move (maximal deformation: " << maximalDeformation
//                << ", lambda: " << m_Lambda << ")" << std::endl;
//      std::cout << "        Cost going from " << currentCost << " to " << trialCost << std::endl;
 //     std::cout << "        Decreasing lambda from " << m_Lambda;

      m_Lambda /= lambdaDecreaseFactor;

//      std::cout << " to " << m_Lambda << std::endl;

      currentCost = trialCost;

      AtlasMeshDeformationLevenbergMarquardt::Pointer  tmp = m_LevenbergMarquardt;
      m_LevenbergMarquardt = m_TrialLevenbergMarquardt;
      m_TrialLevenbergMarquardt = tmp;

      return maximalDeformation;
      }
    else
      {
      // Bad move. Restore previous position
//      std::cout << "Bad trial move (maximal deformation: " << maximalDeformation
//                << ", lambda: " << m_Lambda << ")" << std::endl;
//      std::cout << "        Cost going from " << currentCost << " to " << trialCost << std::endl;

      m_Mesh->SetPoints( currentPosition );

#if 0    
      if ( m_LevenbergMarquardt->GetPeerToCopyHessianFrom() )
        {
        //
        std::cout << "We were using a lazy Hessian - let's compute it properly and try with same lambda" << std::endl;
        m_LevenbergMarquardt->SetPeerToCopyHessianFrom( 0 );
        m_LevenbergMarquardt->Rasterize( m_Mesh );
        }
      else
        {
#endif
        // Increase lambda and try again - unless the proposed deformation was already so small you want to stop altogether
        if ( maximalDeformation > 0.01 /* 0.01 */ )
          {
//          std::cout << "        Increasing lambda from " << m_Lambda;

          m_Lambda *= lambdaIncreaseFactor;

//          std::cout << " to " << m_Lambda << " and trying again." << std::endl;
          }
        else
          {
 //         std::cout << "        proposed deformation is already too small; bailing out" << std::endl;
        
          return 0;
          }
#if 0
        }
#endif
      } // End test if good or bad trial move


    } // End trying loop until we find a good step or we're converged
 

}





} // end namespace kvl

