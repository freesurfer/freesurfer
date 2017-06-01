#include "kvlAtlasMeshDeformationFixedStepGradientDescentOptimizer.h"


namespace kvl
{


//
//
//
AtlasMeshDeformationFixedStepGradientDescentOptimizer
::AtlasMeshDeformationFixedStepGradientDescentOptimizer()
{
  m_StepSize = 1.0;
  m_LineSearchStopCriterion = 0.00005;
}


//
//
//
AtlasMeshDeformationFixedStepGradientDescentOptimizer
::~AtlasMeshDeformationFixedStepGradientDescentOptimizer()
{
}




  
//
//
//
double
AtlasMeshDeformationFixedStepGradientDescentOptimizer
::FindAndOptimizeNewSearchDirection()
{

  // Compute the largest gradient magnitude in any point
  const double  maximumGradientMagnitude = this->ComputeMaximalDeformation( m_Gradient );
  if ( maximumGradientMagnitude == 0 )
    {
    return 0.0;  
    }
    
  // Try to add the scaled gradient to the current position to obtain the trial position
  const double  alpha = -( m_StepSize / maximumGradientMagnitude );
  double  maximalDeformation = 0.0;
  AtlasMesh::PointsContainer::Pointer  trialPosition = 0;
  this->AddDeformation( m_Position, alpha, m_Gradient, trialPosition, maximalDeformation );
  if ( m_Verbose )
    {
    std::cout << "Using trialPosition with maximalDeformation: " << maximalDeformation << std::endl;  
    }
  
  // Try out this new position
  AtlasPositionGradientContainerType::Pointer  trialGradient = 0;
  double  trialCost = 0.0;
  this->GetCostAndGradient( trialPosition, trialCost, trialGradient );
  if ( m_Verbose )
    {
    std::cout << "m_Cost: " << m_Cost << std::endl;  
    std::cout << "trialCost: " << trialCost << std::endl;  
    }
  if ( ( trialCost > m_Cost ) ||
       ( ( fabsf( m_Cost - trialCost ) / fabsf( trialCost ) ) < m_LineSearchStopCriterion ) )
    {
    // Bad or insufficiently good step -- give up
    return 0.0;
    }
    
  // Successful move
  m_Position = trialPosition;
  m_Cost = trialCost;
  m_Gradient = trialGradient;
        
  return maximalDeformation;        
        
}        



} // end namespace kvl

