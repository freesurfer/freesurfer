#include "kvlAtlasMeshDeformationOptimizer.h"


namespace kvl
{


//
//
//
AtlasMeshDeformationOptimizer
::AtlasMeshDeformationOptimizer()
{
  m_Mesh = 0;

  m_IterationNumber = 0;
  m_MaximumNumberOfIterations = 1000 /* itk::NumericTraits< unsigned int >::max() */;
  m_IterationEventResolution = 10;
  m_Initialized = false;
  m_Verbose = false;
  m_Calculator = 0;
  m_MaximalDeformationStopCriterion = 0.05;
  
  m_Cost = 0;
  m_Position = 0;
  m_Gradient = 0;
  
  m_LineSearchMaximalDeformationLimit = 50.0; // Measured in voxels
  m_LineSearchMaximalDeformationIntervalStopCriterion = 0.05; // Measured in voxels
  
}


//
//
//
AtlasMeshDeformationOptimizer
::~AtlasMeshDeformationOptimizer()
{
}


//
//
//
void
AtlasMeshDeformationOptimizer
::GetCostAndGradient( const AtlasMesh::PointsContainer* position, 
                      double& cost, 
                      AtlasPositionGradientContainerType::Pointer& gradient )
{
  
  // Set up a mesh
  AtlasMesh::Pointer  mesh = AtlasMesh::New();
  mesh->SetPoints( const_cast< AtlasMesh::PointsContainer* >( position ) );
  mesh->SetCells( m_Mesh->GetCells() );
  mesh->SetPointData( m_Mesh->GetPointData() );
  mesh->SetCellData( m_Mesh->GetCellData() );
    
  // Rasterize mesh
  m_Calculator->Rasterize( mesh );
  
  // Retrieve results
  cost = m_Calculator->GetMinLogLikelihoodTimesPrior();
  gradient = m_Calculator->GetPositionGradient();

  
}
  
  
  
//
//
//
bool
AtlasMeshDeformationOptimizer
::Go()
{
  bool hasMoved = false;
  while ( true )
    {
    const double  maximalDeformation = this->Step();
    if ( maximalDeformation > 0.0 )
      {
      hasMoved = true;  
      }
    else
      {
      break;
      }  
    }  
  
  return hasMoved;  
}



//
//
//
double
AtlasMeshDeformationOptimizer
::Step()
{

  // Test if we're running out of time
  if ( m_IterationNumber >= m_MaximumNumberOfIterations )
    {
    this->InvokeEvent( DeformationEndEvent() );
    return 0.0;  
    }

  // Make sure m_Position and m_Gradient are available
  if ( !m_Initialized )
    {
    this->Initialize();  
    }
    
    
  // Try to make one successful move
  double maximalDeformation = 0.0;
  try
    {
    maximalDeformation = this->FindAndOptimizeNewSearchDirection();
    }
  catch( itk::ExceptionObject& err )
    {
    std::cout << "Something wrong with finding and/or optimizing the search direction -- giving up" << std::endl; 
    std::cout << "(" << err << ")" << std::endl;  
    this->InvokeEvent( DeformationEndEvent() );
    return 0.0;  
    }
    
  m_Mesh->SetPoints( m_Position ); 
  //std::cout << "m_Mesh->GetPoints()->Begin().Value(): " << m_Mesh->GetPoints()->Begin().Value() << std::endl;
  
  if ( !( m_IterationNumber % m_IterationEventResolution ) )
    {
    this->InvokeEvent( DeformationIterationEvent() );
    }
    
    
  // Test if we need to stop
  if ( maximalDeformation <= m_MaximalDeformationStopCriterion )
    {
    if ( ( m_IterationNumber > 0 ) || ( maximalDeformation == 0.0 ) )
      {
      if ( m_Verbose )
        {
        std::cout << "Optimizer: maximalDeformation is too small; stopping" << std::endl;
        }
      this->InvokeEvent( DeformationEndEvent() );
      return 0.0;  
      }
    }
    
  //
  m_IterationNumber++;
  return maximalDeformation;
    
}


  
//
//
//
double
AtlasMeshDeformationOptimizer
::ComputeMaximalDeformation( const AtlasPositionGradientContainerType* deformation ) const
{
  
  // Compute the largest deformation magnitude in any point
  double  maximalDeformation = 0.0;
  for ( AtlasPositionGradientContainerType::ConstIterator  it = deformation->Begin();
        it != deformation->End(); ++it )
    {     
    const double  magnitude = it.Value().GetNorm(); 
    if ( magnitude > maximalDeformation )
      {
      maximalDeformation = magnitude;
      }
    }

  return maximalDeformation;
}
  
  
  
//
//
//
void
AtlasMeshDeformationOptimizer
::AddDeformation( const AtlasMesh::PointsContainer* position, 
                  double alpha,
                  const AtlasPositionGradientContainerType* deformationDirection,                    
                  AtlasMesh::PointsContainer::Pointer&  newPosition,
                  double&  maximalDeformation ) const
{

  maximalDeformation = 0.0;

  newPosition = AtlasMesh::PointsContainer::New();
  AtlasMesh::PointsContainer::ConstIterator  posIt = position->Begin();
  AtlasPositionGradientContainerType::ConstIterator  defIt = deformationDirection->Begin();
  for ( ; posIt != position->End(); ++posIt, ++defIt )
    {
    const AtlasPositionGradientType  newStep  = alpha * defIt.Value();

    if ( newStep.GetNorm() > maximalDeformation )
      {
      maximalDeformation = newStep.GetNorm();
      }

    newPosition->InsertElement( posIt.Index(), posIt.Value() + newStep );
  
    }
    
}                  
  

  
//
//
//
double
AtlasMeshDeformationOptimizer
::ComputeInnerProduct( const AtlasPositionGradientContainerType* deformation1,
                       const AtlasPositionGradientContainerType* deformation2 ) const
{
  
  double innerProduct = 0.0;
  AtlasPositionGradientContainerType::ConstIterator  defIt1 = deformation1->Begin();
  AtlasPositionGradientContainerType::ConstIterator  defIt2 = deformation2->Begin();
  for ( ; defIt1 != deformation1->End(); ++defIt1, ++defIt2 )
    {
    for ( int i = 0; i < 3; i++ )
      {
      innerProduct += defIt1.Value()[i] * defIt2.Value()[i];
      }
    }  
  
  return innerProduct;
}                       
  
  
//
//
//
AtlasPositionGradientContainerType::Pointer  
AtlasMeshDeformationOptimizer
::LinearlyCombineDeformations( const AtlasPositionGradientContainerType* deformation1,
                               double beta1,
                               const AtlasPositionGradientContainerType* deformation2,
                               double beta2 ) const
{
  
  AtlasPositionGradientContainerType::Pointer  newDeformation = AtlasPositionGradientContainerType::New();
  AtlasPositionGradientContainerType::ConstIterator  defIt1 = deformation1->Begin();
  AtlasPositionGradientContainerType::ConstIterator  defIt2 = deformation2->Begin();
  for ( ; defIt1 != deformation1->End(); ++defIt1, ++defIt2 )
    {
    newDeformation->InsertElement( defIt1.Index(), beta1 * defIt1.Value() + beta2 * defIt2.Value() );
    }
  
  return newDeformation;
}                                


  
//
//
//
AtlasPositionGradientContainerType::Pointer  
AtlasMeshDeformationOptimizer
::ScaleDeformation( const AtlasPositionGradientContainerType* deformation,
                    double beta ) const
{
  
  AtlasPositionGradientContainerType::Pointer  newDeformation = AtlasPositionGradientContainerType::New();
  AtlasPositionGradientContainerType::ConstIterator  defIt = deformation->Begin();
  for ( ; defIt != deformation->End(); ++defIt )
    {
    newDeformation->InsertElement( defIt.Index(), beta * defIt.Value() );
    }
  
  return newDeformation;
  
}                    



//
//
//
void 
AtlasMeshDeformationOptimizer
::Initialize()
{
    
  if ( !m_Calculator )
    {
    itkExceptionMacro( << "Cost and gradient calculator missing!" );
    }

  //
  m_Position = m_Mesh->GetPoints();
  this->GetCostAndGradient( m_Position, m_Cost, m_Gradient );
  
  m_Initialized = true;
  
  this->InvokeEvent( DeformationStartEvent() );
  
}




//
//
//
void 
AtlasMeshDeformationOptimizer
::DoLineSearch( const AtlasMesh::PointsContainer*  startPosition, 
                double  startCost,
                const AtlasPositionGradientContainerType*  startGradient,                    
                const AtlasPositionGradientContainerType*  searchDirection,                    
                double  startAlpha,
                double  c1,
                double  c2,
                AtlasMesh::PointsContainer::Pointer&  newPosition,
                double&  newCost,
                AtlasPositionGradientContainerType::Pointer& newGradient,
                double&  alphaUsed )
{
  
  //
  // Make sure this function can be called with startPosition and newPosition
  // pointing to the same object: internally newPosition will become a new smart pointer 
  // (New()), decreasing the reference count of startPosition if newPosition points to
  // it, potentially setting it to zero which then destructs startPosition before it's
  // being used.
  AtlasMesh::PointsContainer::ConstPointer  tmp = startPosition;
  
  // 
  // Part I: Bracket minimum
  //
  const double  initialAlpha = 0.0;
  const double  initialCost = startCost;
  AtlasPositionGradientContainerType::ConstPointer  initialGradient = startGradient;
  const double  initialDirectionalDerivative = 
    this->ComputeInnerProduct( startGradient, searchDirection ); // gradient' * p;
  if ( initialDirectionalDerivative >= 0 )
    {
    itkExceptionMacro( << "search direction is not a descent direction!" )
    }


  //
  double  previousAlpha = initialAlpha;
  double  previousCost = initialCost;
  AtlasPositionGradientContainerType::ConstPointer  previousGradient = initialGradient;
  double  previousDirectionalDerivative = initialDirectionalDerivative;

  //
  const double  maximalDeformationOfSearchDirection = this->ComputeMaximalDeformation( searchDirection );
  const double  alphaMax = m_LineSearchMaximalDeformationLimit / maximalDeformationOfSearchDirection;
  double  alpha = startAlpha;
  if ( alpha > alphaMax )
    {
    alpha = alphaMax; // min( alpha1, alphaMax );
    }
  
  bool  bracketingFailed = true;
  const int  maximumOfBracketingIterations = 10;
  double  lowAlpha = 0.0;
  double  lowCost = 0.0;
  AtlasPositionGradientContainerType::ConstPointer  lowGradient = 0;
  double  lowDirectionalDerivative = 0.0;
  double  highAlpha = 0.0;
  double  highCost = 0.0;
  AtlasPositionGradientContainerType::ConstPointer  highGradient = 0;
  double  highDirectionalDerivative = 0.0;
      
  for ( int bracketingIterationNumber = 0; 
        bracketingIterationNumber < maximumOfBracketingIterations; 
        bracketingIterationNumber++ )
    {     
    // Evaluate current alpha: [ cost gradient ] = tryFunction( x + alpha * p );
    AtlasMesh::PointsContainer::Pointer  position = 0;
    double  maximalDeformation = 0.0;
    this->AddDeformation( startPosition, alpha, searchDirection, 
                          position, maximalDeformation );
    double  cost;
    AtlasPositionGradientContainerType::Pointer  gradient = 0;
    this->GetCostAndGradient( position, cost, gradient );
    const double  directionalDerivative = this->ComputeInnerProduct( gradient, searchDirection ); // gradient' * p


    //
    // Check what sort of solution we have: 
    //   * aweful (good news because then we can stop bracketing and move on to zooming); 
    //   * excellent (good news because then we have found a minimum even without zooming);
    //   * OK but with a positive directional derivative (good news because then we can stop 
    //     bracketing and move on to zooming);
    //   * simply OK (bad news because then our bracket is not yet large enough -- i.e., minimum 
    //     might lie farther away then what we've tested so far -- so we need to keep bracketing)
    //
    if ( (  cost > ( initialCost + c1 * alpha * initialDirectionalDerivative ) ) ||
         ( ( cost >= previousCost ) && ( bracketingIterationNumber > 0 ) ) )
      {
      // Found a large enough interval that we get an aweful function value
      if ( m_Verbose )
        {
        std::cout << "[BRACKETING] Hit a really bad solution -- ready for zooming: "
                  << maximalDeformation << std::endl;
        }          
      
      lowAlpha = previousAlpha;
      lowCost = previousCost;
      lowGradient = previousGradient;
      lowDirectionalDerivative = previousDirectionalDerivative;
      
      highAlpha = alpha;
      highCost = cost;
      highGradient = gradient;
      highDirectionalDerivative = directionalDerivative;
      bracketingFailed = false;
      break;
      }

    if ( fabsf( directionalDerivative ) <= ( -c2 * initialDirectionalDerivative ) )
      {
      // Found an excellent solution that we can simply return -- no need for zooming
      newPosition = position;
      newCost = cost;
      newGradient = gradient;
      alphaUsed = alpha;
      if ( m_Verbose )
        {
        std::cout << "[BRACKETING] Found a really good solution -- no need for zooming: "
                  << maximalDeformation << std::endl;
        std::cout << "-----------------------" << std::endl;
        }
      return;
      }

    if ( directionalDerivative > 0 )
      {
      // Neither excellent nor aweful soution; however directional deriviate is uphill here so 
      // we can start zooming 
      if ( m_Verbose )
        {
        std::cout << "[BRACKETING] Found an OK solution that is sufficiently bad for zooming: "
                  << maximalDeformation << std::endl;
        }          
      lowAlpha = alpha;
      lowCost = cost;
      lowGradient = gradient;
      lowDirectionalDerivative = directionalDerivative;
      
      highAlpha = previousAlpha;
      highCost = previousCost;
      highGradient = previousGradient;
      highDirectionalDerivative = previousDirectionalDerivative;
      bracketingFailed = false;
      break;
      }

      
    // If we're here, it means our solution isn't sufficiently bad yet for zooming - 
    // expand in an effort to get something worse
    if ( alpha >= alphaMax )
      {
      if ( m_Verbose )
        {
        std::cout << "[BRACKETING] Would need to expand but can't anymore" << std::endl;
        }
      previousAlpha = alpha;
      previousCost = cost;
      previousGradient = gradient;
      previousDirectionalDerivative = directionalDerivative;
      break;
      }
    
    if ( m_Verbose )
      {
      std::cout << "[BRACKETING] Expanding from: " << maximalDeformation << std::endl;
      }
      
    // If a fitted cubic polynomial has a minimum ("cubic interpolation") use it,
    // but make sure to expand sufficiently to reach alphaMax within maximumOfBracketingIterations
    // iterations (and never beyond alphaMax
#if 0
    const double  allowableRangeMin = exp( log( startAlpha ) + 
                                           ( log( alphaMax ) - log( startAlpha ) ) /
                                           maximumOfBracketingIterations * 
                                           ( bracketingIterationNumber + 1 ) ); 
#else
    const double  allowableRangeMin = exp( log( alpha ) + 
                                           ( log( alphaMax ) - log( alpha ) ) /
                                           ( maximumOfBracketingIterations - bracketingIterationNumber ) ); 
#endif
    const double  allowableRangeMax = alphaMax;
    const double  d1 = previousDirectionalDerivative + directionalDerivative - 
                       3 * ( previousCost - cost ) / ( previousAlpha - alpha );
    const double  d2Square = pow( d1, 2 ) - previousDirectionalDerivative * directionalDerivative;
    double  newAlpha = 0.0;
    if ( d2Square > 0 )
      {
      // Interpolation seems to be possible -- use it unless it's too big or too small
      const double  d2 = sqrt( d2Square );
      newAlpha = alpha - ( alpha - previousAlpha ) * 
                 ( directionalDerivative + d2 - d1 ) /
                 ( directionalDerivative - previousDirectionalDerivative + 2 * d2 );
      if ( m_Verbose )
        {
        std::cout << "   interpolation suggestion: " 
                  << newAlpha * maximalDeformationOfSearchDirection << std::endl;
        }          
      
      if ( newAlpha < allowableRangeMin )
        {
        newAlpha = allowableRangeMin;
        if ( m_Verbose )
          {
          std::cout << "   clipped from left: " 
                    << newAlpha * maximalDeformationOfSearchDirection << std::endl;
          }           
        }
      if ( newAlpha > allowableRangeMax )
        {
        newAlpha = allowableRangeMax;
        if ( m_Verbose )
          {
          std::cout << "   clipped from right: " 
                    << newAlpha * maximalDeformationOfSearchDirection << std::endl;
          }          
        }
      }  
    else
      {
      newAlpha = allowableRangeMin;
      if ( m_Verbose )
        {
        std::cout << "   forced to guess work: " 
                  << newAlpha * maximalDeformationOfSearchDirection << std::endl;
        }          
      }

      
    previousAlpha = alpha;
    previousCost = cost;
    previousGradient = gradient;
    previousDirectionalDerivative = directionalDerivative;
    
    alpha = newAlpha;
    
    } // End loop over bracketing iterations

  
  // 
  if ( bracketingFailed )
    {
    // Bracketing failed -- never found an end point that was either bad enough
    // or was magically a good solution itself. 
    if ( m_Verbose )
      {
      std::cout << "Line search couldn't find a suitable bracketing interval" << std::endl;
      }
    
    double  dummy = 0.0;
    this->AddDeformation( startPosition, previousAlpha, searchDirection, 
                          newPosition, dummy ); // xStar = x + previousAlpha * p;
    newCost = previousCost;
    newGradient = const_cast< AtlasPositionGradientContainerType* >( previousGradient.GetPointer() );
    alphaUsed = previousAlpha;
    if ( m_Verbose )
      {
      std::cout << "-----------------------" << std::endl;
      }
    return;
    }

  
  
  //
  // Part II: now that we're here, we have an interval so we can "zoom" 
  //
  if ( m_Verbose )
    {
    std::cout << "[ZOOMING] Zooming with low=" 
              << lowAlpha * maximalDeformationOfSearchDirection 
              << " and high=" 
              << highAlpha * maximalDeformationOfSearchDirection 
              << std::endl;
    }          

  while ( true )
    {
    // If possibly use minimum of cubic polynomial approximation of function; otherwise bisect
    double  leftAlpha = 0.0;
    double  leftCost = 0.0;
    double  leftDirectionalDerivative = 0.0;
  
    double  rightAlpha = 0.0;
    double  rightCost = 0.0;
    double  rightDirectionalDerivative = 0.0;
      
    if ( highAlpha > lowAlpha )
      {
      rightAlpha = highAlpha;
      rightCost = highCost;
      rightDirectionalDerivative = highDirectionalDerivative;
      
      leftAlpha = lowAlpha;
      leftCost = lowCost;
      leftDirectionalDerivative = lowDirectionalDerivative;
      }
    else
      {
      leftAlpha = highAlpha;
      leftCost = highCost;
      leftDirectionalDerivative = highDirectionalDerivative;
      
      rightAlpha = lowAlpha;
      rightCost = lowCost;
      rightDirectionalDerivative = lowDirectionalDerivative;
      }
    
    const double  d1 = leftDirectionalDerivative + rightDirectionalDerivative - 
                       3 * ( leftCost - rightCost ) / ( leftAlpha - rightAlpha );
    const double  d2Square = pow( d1, 2 ) - leftDirectionalDerivative * rightDirectionalDerivative;
    if ( ( d2Square < 0 ) || ( highCost == itk::NumericTraits< double >::max() ) )
      {
      // 
      alpha = ( leftAlpha + rightAlpha ) / 2;
      if ( m_Verbose )
        {
        std::cout << "   bisected: " << alpha * maximalDeformationOfSearchDirection << std::endl;
        // std::cout << "         [ bisected because lowCost=" << lowCost 
        //           << " and highCost=" << highCost << "]" << std::endl;
        }
      }
    else
      {
      const double  d2 = sqrt( d2Square );
      alpha = rightAlpha - ( rightAlpha - leftAlpha ) * 
              ( rightDirectionalDerivative + d2 - d1 ) / 
              ( rightDirectionalDerivative -leftDirectionalDerivative + 2 * d2 );
      if ( m_Verbose )
        {
        std::cout << "   interpolated: " << alpha * maximalDeformationOfSearchDirection << std::endl;
        }
      }

      
    // Make sure alpha is clipped to area well within [ lowAlpha, highAlpha ] range
    double  acceptableRangeMin = 0.0;
    double  acceptableRangeMax = 0.0;
    if ( highAlpha > lowAlpha )
      {
      acceptableRangeMin = lowAlpha + 0.1 * ( highAlpha - lowAlpha );
      acceptableRangeMax = highAlpha - 0.1 * ( highAlpha - lowAlpha );     
      }
    else
      {
      acceptableRangeMin = highAlpha + 0.1 * ( lowAlpha - highAlpha );
      acceptableRangeMax = lowAlpha - 0.1 * ( lowAlpha - highAlpha );
      }

#if 0      
    if ( alpha < acceptableRangeMin )
      {
      alpha = ( acceptableRangeMin + acceptableRangeMax ) / 2;
      std::cout << "    Too small for current range [" 
                << leftAlpha * maximalDeformationOfSearchDirection 
                << " " 
                << rightAlpha * maximalDeformationOfSearchDirection 
                << "]; using middle instead: "
                << alpha * maximalDeformationOfSearchDirection << std::endl;
      }
    else if ( alpha > acceptableRangeMax )
      {
      alpha = ( acceptableRangeMin + acceptableRangeMax ) / 2;
      std::cout << "    Too big for current range [" 
                << leftAlpha * maximalDeformationOfSearchDirection 
                << " " 
                << rightAlpha * maximalDeformationOfSearchDirection 
                << "]; using middle instead: "
                << alpha * maximalDeformationOfSearchDirection << std::endl;
      }
#else
    if ( alpha < acceptableRangeMin )
      {
      alpha = acceptableRangeMin;
      if ( m_Verbose )
        {
        std::cout << "    Too small for current range [" 
                  << leftAlpha * maximalDeformationOfSearchDirection 
                  << " " 
                  << rightAlpha * maximalDeformationOfSearchDirection 
                  << "]; using something forcibly bigger instead: "
                  << alpha * maximalDeformationOfSearchDirection << std::endl;
        }          
      }
    else if ( alpha > acceptableRangeMax )
      {
      alpha = acceptableRangeMax;
      if ( m_Verbose )
        {
        std::cout << "    Too big for current range [" 
                  << leftAlpha * maximalDeformationOfSearchDirection 
                  << " " 
                  << rightAlpha * maximalDeformationOfSearchDirection 
                  << "]; using something forcibly smaller instead: "
                  << alpha * maximalDeformationOfSearchDirection << std::endl;
        }            
      }
#endif
    
    // Evaluate cost function: [ cost gradient ] = tryFunction( x + alpha * p );
    AtlasMesh::PointsContainer::Pointer  position = 0;
    double  maximalDeformation = 0.0;
    this->AddDeformation( startPosition, alpha, searchDirection, 
                          position, maximalDeformation );
    double  cost;
    AtlasPositionGradientContainerType::Pointer  gradient = 0;
    this->GetCostAndGradient( position, cost, gradient );
    const double  directionalDerivative = this->ComputeInnerProduct( gradient, searchDirection ); // gradient' * p

    if ( ( cost > ( initialCost + c1 * alpha * initialDirectionalDerivative ) ) || 
         ( cost >= lowCost ) )
      {
      highAlpha = alpha;
      highCost = cost;
      highGradient = gradient;
      highDirectionalDerivative = directionalDerivative;
      }
    else
      {
      if ( abs( directionalDerivative ) <= ( -c2 * initialDirectionalDerivative ) )
        {
        newPosition = position;
        newCost = cost;
        newGradient = gradient;
        alphaUsed = alpha;
        if ( m_Verbose )
          {
          std::cout << "[ZOOMING] Found a good solution: " << maximalDeformation << std::endl;
          std::cout << "-----------------------" << std::endl;
          }
        return;
        }
      
      if ( directionalDerivative * ( highAlpha - lowAlpha ) >= 0 )
        {
        highAlpha = lowAlpha;
        highCost = lowCost;
        highGradient = lowGradient;
        highDirectionalDerivative = lowDirectionalDerivative;
        }
      
      //
      lowAlpha = alpha;
      lowCost = cost;
      lowGradient = gradient;
      lowDirectionalDerivative = directionalDerivative;
      
      }

    
    // Check size
    if ( ( fabsf( highAlpha - lowAlpha ) * maximalDeformationOfSearchDirection )
         < m_LineSearchMaximalDeformationIntervalStopCriterion )
      {
      //std::cout << "!!!!!!!!!!!!!!!" << startPosition->Begin().Value() << std::endl;
      //std::cout << "!!!!!!!!!!!!!!!" << searchDirection->Begin().Value() << std::endl;
      //std::cout << "!!!!!!!!!!!!!!!" << lowAlpha << std::endl;

      double  dummy = 0.0;
      this->AddDeformation( startPosition, lowAlpha, searchDirection, 
                            newPosition, dummy ); // xStar = x + lowAlpha * p;
      
      //std::cout << "!!!!!!!!!!!!!!!" << newPosition->Begin().Value() << std::endl;
      
      newCost = lowCost;
      newGradient = const_cast< AtlasPositionGradientContainerType* >( lowGradient.GetPointer() );
      alphaUsed = lowAlpha;
      if ( m_Verbose )
        {
        std::cout << "[ZOOMING] Zooming interval size too low; returning best: "
                  << lowAlpha * maximalDeformationOfSearchDirection << std::endl;
        std::cout << "-----------------------" << std::endl;        
        }
      return;
      }
    
    } // End zoom iterations

  
}


                       
} // end namespace kvl

