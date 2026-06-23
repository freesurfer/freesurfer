#include "kvlAtlasMeshDeformationGradientDescentOptimizer.h"


namespace kvl
{


//
//
//
AtlasMeshDeformationGradientDescentOptimizer
::AtlasMeshDeformationGradientDescentOptimizer()
{

  m_OldCost = 0;
  m_OldGradient = 0;
  m_OldSearchDirection = 0;
  m_AlphaUsedLastTime = 0.0;
  
  m_StartDistance = 1.0; // Measured in voxels

}


//
//
//
AtlasMeshDeformationGradientDescentOptimizer
::~AtlasMeshDeformationGradientDescentOptimizer()
{
}



//
//
//
void AtlasMeshDeformationGradientDescentOptimizer
::Initialize()
{

  m_OldCost = 0;
  m_OldGradient = 0;
  m_OldSearchDirection = 0;
  
  Superclass::Initialize();
  
}

  
  
//
//
//
double
AtlasMeshDeformationGradientDescentOptimizer
::FindAndOptimizeNewSearchDirection()
{

  // 
  // Part I: Decide on a new search direction
  //
  AtlasPositionGradientContainerType::Pointer  searchDirection 
              = this->ScaleDeformation( m_Gradient, -1.0 );   // p = -gradient;


  //
  // PartII: Make an educated guess of the appropriate step size
  //
  double  startAlpha = 0.0;
  if ( this->GetIterationNumber() == 0 )
    {
    startAlpha = m_StartDistance / this->ComputeMaximalDeformation( searchDirection );
    }
  else
    {
    if ( false )
      {
      // Do something similar to previous iteration
      //alpha1 = alphaUsed * ( gradientOld' * pOld ) / ( gradient' * p );
      startAlpha = m_AlphaUsedLastTime * 
                   this->ComputeInnerProduct( m_OldGradient, m_OldSearchDirection ) /
                   this->ComputeInnerProduct( m_Gradient, searchDirection );
      }
    else if ( false )
      {
      // Use optimum of quadratic approximation
      //alpha1 = 2 * ( cost - costOld ) / ( gradient' * p );
      startAlpha = 2 * ( m_Cost - m_OldCost ) / 
                   this->ComputeInnerProduct( m_Gradient, searchDirection );
      }
    else
      {
      // Use domain knowledge: try same step size in terms of number of voxels
      // as last time
      //alpha1 = alphaUsed * max( abs( pOld ) ) / max( abs( p ) ); 
      startAlpha = m_AlphaUsedLastTime * 
                   this->ComputeMaximalDeformation( m_OldSearchDirection ) / 
                   this->ComputeMaximalDeformation( searchDirection );
      }
    } // End test if starting/restarting
  

  
  //
  // Part III: Line Search
  // 
  const double  c1 = 1e-4;
  const double  c2 = 0.9;
  m_OldCost = m_Cost;
  m_OldGradient = m_Gradient;
  m_OldSearchDirection = searchDirection;
  // [ x, cost, gradient, alphaUsed ] = tryLineSearch( x, cost, gradient, p, alpha1, c1, c2 );
  double  alphaUsed = 0.0;
  this->DoLineSearch( m_Position, 
                      m_Cost,
                      m_Gradient,                    
                      searchDirection,                    
                      startAlpha,
                      c1,
                      c2,
                      m_Position,
                      m_Cost,
                      m_Gradient,
                      alphaUsed );      
  
  //std::cout << "m_Cost: " << m_Cost << std::endl;


  // Some book keeping
  const double  maximalDeformation = alphaUsed * this->ComputeMaximalDeformation( searchDirection );
  m_AlphaUsedLastTime = alphaUsed;
  
  return maximalDeformation;
}







} // end namespace kvl

