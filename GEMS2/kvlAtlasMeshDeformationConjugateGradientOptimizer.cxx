#include "kvlAtlasMeshDeformationConjugateGradientOptimizer.h"


namespace kvl
{


//
//
//
AtlasMeshDeformationConjugateGradientOptimizer
::AtlasMeshDeformationConjugateGradientOptimizer()
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
AtlasMeshDeformationConjugateGradientOptimizer
::~AtlasMeshDeformationConjugateGradientOptimizer()
{
}



//
//
//
void AtlasMeshDeformationConjugateGradientOptimizer
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
AtlasMeshDeformationConjugateGradientOptimizer
::FindAndOptimizeNewSearchDirection()
{

  // 
  // Part I: Decide on a new search direction
  //
  AtlasPositionGradientContainerType::Pointer  searchDirection = 0;
  bool  startingOrRestarting = false;
  if ( this->GetIterationNumber() == 0 )
    {
    startingOrRestarting = true;
    }
  else if ( ( this->ComputeInnerProduct( m_Gradient, m_OldGradient ) / 
              this->ComputeInnerProduct( m_Gradient, m_Gradient ) ) 
             >= 0.1 ) // ( ( gradient' * gradientOld ) / ( gradient' * gradient ) >= 0.1 )
    {
    // Somehow we didn't really get rid of the previous direction; we need a restart
    if ( m_Verbose )
      {
      std::cout << "Somehow didn't get rid of previous direction; restarting" << std::endl;
      }
    startingOrRestarting = true;
    }

  if ( !startingOrRestarting )
    {
    // Normal regime
    double  beta = 0.0;
    if ( false )
      {
      // Polak-Ribiere
      // beta = gradient' * ( gradient - gradientOld ) / ( gradientOld' * gradientOld ); 
      // beta = max( beta, 0 );
      beta = this->ComputeInnerProduct( m_Gradient,
                                        this->LinearlyCombineDeformations( m_Gradient, 1.0, 
                                                                           m_OldGradient, -1.0 ) ) /
             this->ComputeInnerProduct( m_OldGradient, m_OldGradient ); 

      if ( beta < 0 )
        {
        beta = 0.0;
        }
        
      }
    else
      {
      // Hestenes-Stiefel
      // beta = gradient' * ( gradient - gradientOld ) / ( ( gradient - gradientOld )' * pOld ); 
      AtlasPositionGradientContainerType::Pointer  tmp =  
              this->LinearlyCombineDeformations( m_Gradient, 1.0, 
                                                 m_OldGradient, -1.0 );
      beta = this->ComputeInnerProduct( m_Gradient, tmp ) / 
             this->ComputeInnerProduct( tmp, m_OldSearchDirection );       
      }
      
      
    // p = -gradient + beta * pOld;
    searchDirection = this->LinearlyCombineDeformations( m_Gradient, -1.0, m_OldSearchDirection, beta );
    
    
    if ( this->ComputeInnerProduct( searchDirection, m_Gradient ) > 0 ) // ( p' * gradient > 0 )
      {
      // Senseless direction proposed - restart anyway
      if ( m_Verbose )
        {
        std::cout << "Proposed direction is not a descent direction; restarting" << std::endl;
        }
      startingOrRestarting = true;
      }

    } // End testing if we need to start/restart

  if ( startingOrRestarting )
    {
    if ( m_Verbose )
      {
      std::cout << "(Re)starting conjugate gradient" << std::endl;
      }

    // p = -gradient;
    searchDirection = this->ScaleDeformation( m_Gradient, -1.0 ); 
    }

  
  
  
  //
  // PartII: Make an educated guess of the appropriate step size
  //
  double  startAlpha = 0.0;
  if ( startingOrRestarting )
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
  const double  c2 = 0.1; // Accurate line search required from conjugate gradient
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

