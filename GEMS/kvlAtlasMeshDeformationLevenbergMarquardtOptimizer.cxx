#include "kvlAtlasMeshDeformationLevenbergMarquardtOptimizer.h"


namespace kvl
{


//
//
//
AtlasMeshDeformationLevenbergMarquardtOptimizer
::AtlasMeshDeformationLevenbergMarquardtOptimizer()
{
  m_Image = 0;
  m_ProbabilityImage = 0;
  m_Mesh = 0;
  m_MeshToImageTransform = 0;

  m_Initialized = false;
  m_UseProbabilityImage = false;

  m_LevenbergMarquardt = AtlasMeshDeformationLevenbergMarquardt::New();
  m_TrialLevenbergMarquardt = AtlasMeshDeformationLevenbergMarquardt::New();

  m_Lambda = 1.0f;
  
  m_PositionUpdatingIterationNumber = 0;
  m_PositionUpdatingMaximumNumberOfIterations = 11 /* 2000 */ /* itk::NumericTraits< unsigned int >::max() */;
  m_PositionUpdatingIterationEventResolution = 10;
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
  
  // Create an empty, dummy template image needed by anything that rasterizes meshes
  LabelImageType::Pointer  dummyTemplateImage = LabelImageType::New();
  if ( m_UseProbabilityImage )
    {
    dummyTemplateImage->SetRegions( m_ProbabilityImage->GetBufferedRegion() );
    }
  else
    {
    dummyTemplateImage->SetRegions( m_Image->GetBufferedRegion() );
    }
  dummyTemplateImage->Allocate();

  
  // Set up m_LevenbergMarquardt
  m_LevenbergMarquardt->SetLabelImage( dummyTemplateImage );
  m_LevenbergMarquardt->SetUseProbabilityImage( m_UseProbabilityImage );
  m_LevenbergMarquardt->SetImage( m_Image );
  m_LevenbergMarquardt->SetProbabilityImage( m_ProbabilityImage );
  m_LevenbergMarquardt->SetMeans( m_Means );
  m_LevenbergMarquardt->SetVariances( m_Variances );
  m_LevenbergMarquardt->SetMeshToImageTransform( m_MeshToImageTransform );
  
  // Idem for m_TrialLevenbergMarquardt
  m_TrialLevenbergMarquardt->SetLabelImage( dummyTemplateImage );
  m_TrialLevenbergMarquardt->SetUseProbabilityImage( m_UseProbabilityImage );
  m_TrialLevenbergMarquardt->SetImage( m_Image );
  m_TrialLevenbergMarquardt->SetProbabilityImage( m_ProbabilityImage );
  m_TrialLevenbergMarquardt->SetMeans( m_Means );
  m_TrialLevenbergMarquardt->SetVariances( m_Variances );
  m_TrialLevenbergMarquardt->SetMeshToImageTransform( m_MeshToImageTransform );

  // Rasterize m_LevenbergMarquardt
  m_LevenbergMarquardt->Rasterize( m_Mesh );
  
  //m_Lambda = 1.0f;
  m_Initialized = true;
  
}



//
//
//
bool
AtlasMeshDeformationLevenbergMarquardtOptimizer
::Go()
{

  //this->InvokeEvent( PositionUpdatingStartEvent() );

  std::cout << "m_PositionUpdatingMaximumNumberOfIterations: " << m_PositionUpdatingMaximumNumberOfIterations << std::endl;

  
  bool haveMoved = false; // Keep track if we've ever moved or not
  for ( m_PositionUpdatingIterationNumber = 0; 
        m_PositionUpdatingIterationNumber < m_PositionUpdatingMaximumNumberOfIterations;
        m_PositionUpdatingIterationNumber++ )
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
      std::cout << "        maximalDeformation is too small; stopping" << std::endl;
      break;
      }

    if ( !( m_PositionUpdatingIterationNumber % m_PositionUpdatingIterationEventResolution ) )
      {
      //this->InvokeEvent( PositionUpdatingIterationEvent() );
      }

    } // End loop over all positive moves


  //this->InvokeEvent( PositionUpdatingEndEvent() );


  return haveMoved;
 
}



//
//
//
float
AtlasMeshDeformationLevenbergMarquardtOptimizer
::PerformOneSuccessfulStep()
{
  
  const float  lambdaIncreaseFactor = 2.0f; // Factor by which lambda is increased if failure detected
  const float  lambdaDecreaseFactor = 1.1f; // Factor by which lambda is decreased if success detected

  
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
    std::cout << "Calculating Levenberg-Marquardt step with lambda " << m_Lambda << "..." << std::endl;
    AtlasPositionGradientContainerType::Pointer  trialStep = m_LevenbergMarquardt->GetStep( m_Lambda, false );
    std::cout << "...done!" << std::endl;

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

    std::cout << "maximalDeformation: " << maximalDeformation << std::endl;


    
    
    // Evaluate the trial position
    m_Mesh->SetPoints( trialPosition );
    m_TrialLevenbergMarquardt->Rasterize( m_Mesh );
    double  trialCost = m_TrialLevenbergMarquardt->GetMinLogLikelihoodTimesPrior();


    // Do the Right Thing
    std::cout << "=======================" << std::endl;
    if ( trialCost < currentCost )
      {
      // Good move. Keep the move and decrease Marquardt's lambda
      std::cout << "Good trial move (maximal deformation: " << maximalDeformation
                << ", lambda: " << m_Lambda << ")" << std::endl;
      std::cout << "        Cost going from " << currentCost << " to " << trialCost << std::endl;
      std::cout << "        Decreasing lambda from " << m_Lambda;
      m_Lambda /= lambdaDecreaseFactor;
      std::cout << " to " << m_Lambda << std::endl;

      currentCost = trialCost;
      
      AtlasMeshDeformationLevenbergMarquardt::Pointer  tmp = m_LevenbergMarquardt;
      m_LevenbergMarquardt = m_TrialLevenbergMarquardt;
      m_TrialLevenbergMarquardt = tmp;

      return maximalDeformation;
      }
    else
      {
      // Bad move. Restore previous position
      std::cout << "Bad trial move (maximal deformation: " << maximalDeformation
                << ", lambda: " << m_Lambda << ")" << std::endl;
      std::cout << "        Cost going from " << currentCost << " to " << trialCost << std::endl;

      m_Mesh->SetPoints( currentPosition );

      // Increase lambda and try again - unless the proposed deformation was already so small you want to stop altogether
      if ( maximalDeformation > 0.01 /* 0.01 */ )
        {
        std::cout << "        Increasing lambda from " << m_Lambda;
        m_Lambda *= lambdaIncreaseFactor;
        std::cout << " to " << m_Lambda << " and trying again." << std::endl;
        }
      else
        {
        std::cout << "        proposed deformation is already too small; bailing out" << std::endl;
      
        return 0;
        }

      } // End test if good or bad trial move


    } // End trying loop until we find a good step or we're converged

  
  
}





} // end namespace kvl

