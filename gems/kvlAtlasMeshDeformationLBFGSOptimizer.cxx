#include "kvlAtlasMeshDeformationLBFGSOptimizer.h"

#define KVL_ENABLE_TIME_PROBE3 0

#if KVL_ENABLE_TIME_PROBE3
  #include "itkTimeProbe.h"
#endif


namespace kvl
{


//
//
//
AtlasMeshDeformationLBFGSOptimizer
::AtlasMeshDeformationLBFGSOptimizer()
{

  m_OldCost = 0;
  m_OldGradient = 0;
  m_OldSearchDirection = 0;
  m_AlphaUsedLastTime = 0.0;
  
  m_StartDistance = 1.0; // Measured in voxels

  m_MaximumMemoryLength = 12; // Less than 20 recommended
}


//
//
//
AtlasMeshDeformationLBFGSOptimizer
::~AtlasMeshDeformationLBFGSOptimizer()
{
}



//
//
//
void AtlasMeshDeformationLBFGSOptimizer
::WipeMemory()
{
  
  m_OldCost = 0;
  m_OldGradient = 0;
  m_OldSearchDirection = 0;
  
  m_Ss.clear();
  m_Ys.clear();
  m_InverseRhos.clear();
  
}

  
  
//
//
//
double
AtlasMeshDeformationLBFGSOptimizer
::FindAndOptimizeNewSearchDirection()
{

  //
  // BFGS does
  // 
  //    newPosition = oldPosition - alpha * H * gradient
  //
  // where H is an approximation of the inverse Hessian in two steps:
  // 
  // * Part I: Compute the search direction -H*gradient, by first 
  //           updating H based on the estimate we used in the previous
  //           iteration and the change in position and gradient in that
  //           iteration, and then computing -H*gradient
  //
  // * Part II: Do a line search to find the scalar alpha to determine
  //            the correct setp size
  // 
  // We initialize H in the very first iteration as H_0 = gamma * I (diagonal matrix),
  // where gamma is computed so that H_0 * gradient makes a change in position of a
  // physically plausable size (a value alpha=1 is the first guestimate of the line
  // search in the very first iteration, so you're simply doing simple gradient descent
  // of a reasonable step size in the very first iteration).
  //
  // Now L-BFGS does three smart things: 
  //   (1) instead of storing H and gradually updating it in each iteration, it always 
  //       starts from scratch from H_0 and makes the relevant computations on-the-fly 
  //       by memorizing the differences in position and gradient across previous iterations
  //       (up to a fixed, smallish number, forgetting old iterations); 
  //   (2) instead of computing H explicitly and then computing H * gradient, the two
  //       processes are intertwined so you never have to explicitly handle the massive 
  //       full H matrix
  //   (3) sneakily, instead of using always the same H_0 = gamma * I, it uses a different
  //       gamma for each new iteration, determined with some prescribed recipe
  
  
#if KVL_ENABLE_TIME_PROBE3  
  itk::TimeProbe clock;
  clock.Start();
#endif
  
  
  // 
  // Part I: Decide on a new search direction
  //

  // Compute r = H * gradient without ever computing H explicitly    
  double  gamma = 0.0;
  if ( m_OldSearchDirection.GetPointer() == 0 )
    {
    // Make sure that first try of line search will be given 
    // by distance provided by user (first iteration means
    // p = -gradient, and alpha1 of line search is always 1.0 
    // for L-BFGS
    //gamma = initialAlpha1Distance / max( abs( gradient ) ); 
    gamma = m_StartDistance / this->ComputeMaximalDeformation( m_Gradient );
    }
  else
    {
    // Update S and Y in L-BFGS
    // s = x - xOld; % Equivalent to alphaUsed * pOld;
    AtlasPositionGradientContainerType::ConstPointer  s  
                     = this->ScaleDeformation( m_OldSearchDirection, m_AlphaUsedLastTime ).GetPointer();
      
    // y = gradient - gradientOld;
    AtlasPositionGradientContainerType::ConstPointer  y 
                     =  this->LinearlyCombineDeformations( m_Gradient, 1.0, m_OldGradient, -1.0 ).GetPointer();
                     
    // inverseRho = ( s' * y );
    const double  inverseRho = this->ComputeInnerProduct( s, y );
  
    if ( inverseRho > 1e-10 )
      {
      // Add in front 
      m_Ss.insert( m_Ss.begin(), s );
      m_Ys.insert( m_Ys.begin(), y );
      m_InverseRhos.insert( m_InverseRhos.begin(), inverseRho );
    
      if ( m_Ss.size() > m_MaximumMemoryLength )
        {
        // Forget the oldest information (i.e., dump the last element)
        m_Ss.pop_back();
        m_Ys.pop_back();
        m_InverseRhos.pop_back();
        }
  
      // gamma = ( s' * y ) / ( y' * y );
      gamma = this->ComputeInnerProduct( s, y ) / 
              this->ComputeInnerProduct( y, y );
      }
    else
      {
      if ( m_Verbose )
        {
        std::cout << "Skipped L-BFGS history update" << std::endl;
        }
      }
    
    } // End test if first iteration
    
  
  // q = gradient;
  AtlasPositionGradientContainerType::Pointer  q = this->ScaleDeformation( m_Gradient, 1.0 );
  const int  memoryLength = m_Ss.size();
  //std::cout << "memoryLength: " << memoryLength << std::endl;
  
  std::vector< double >  alps( memoryLength, 0.0 );
  for ( int i = 0; i < memoryLength; i++ )
    {
    AtlasPositionGradientContainerType::ConstPointer  s = m_Ss[ i ];
    AtlasPositionGradientContainerType::ConstPointer  y = m_Ys[ i ];
    const double  inverseRho = m_InverseRhos[ i ];
  
    // alp = ( s' * q ) / inverseRho;
    const double  alp = this->ComputeInnerProduct( s, q ) / inverseRho;
  
    // q = q - alp * y;  
    q = LinearlyCombineDeformations( q, 1.0, y, -alp );
  
    alps[ i ] = alp;
    }
    
  // r = gamma * q;
  AtlasPositionGradientContainerType::Pointer  r = this->ScaleDeformation( q, gamma );
  for ( int i = ( memoryLength-1 );  i >=0; i-- )
    {
    AtlasPositionGradientContainerType::ConstPointer  s = m_Ss[ i ];
    AtlasPositionGradientContainerType::ConstPointer  y = m_Ys[ i ];
    const double  inverseRho = m_InverseRhos[ i ];
    const double  alp = alps[ i ];
  
    // bet = ( y' * r ) / inverseRho;
    const double  bet = this->ComputeInnerProduct( y, r ) / inverseRho;

    // r = r + s * ( alp - bet );
    r = this->LinearlyCombineDeformations( r, 1.0, s, alp - bet );
    }
  
  // Direction is -r: p = -r;
  AtlasPositionGradientContainerType::Pointer  searchDirection 
                                                      = this->ScaleDeformation( r, -1.0 );

                                                      
#if KVL_ENABLE_TIME_PROBE3     
  clock.Stop();
  std::cout << "  --- Time taken by determining search direction: " << clock.GetMean() << std::endl;
  clock.Reset();
  clock.Start();
#endif 
                                                
                                                      
                                                      
  //
  // PartII: Make an educated guess of the appropriate step size
  //
  const double  startAlpha = 1.0; // BFGS requires this to proof for convergence properties 
  

  
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

#if KVL_ENABLE_TIME_PROBE3     
  clock.Stop();
  std::cout << "  --- Time taken by line search: " << clock.GetMean() << std::endl;
#endif 
  

  // Some book keeping
  const double  maximalDeformation = alphaUsed * this->ComputeMaximalDeformation( searchDirection );
  m_AlphaUsedLastTime = alphaUsed;
  
  return maximalDeformation;
}







} // end namespace kvl

