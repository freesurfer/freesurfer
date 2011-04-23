#ifndef _itkPowellOptimizer_txx
#define _itkPowellOptimizer_txx

#include "itkPowellOptimizer.h"
#include "itkCommand.h"
#include "itkEventObject.h"
#include "itkExceptionObject.h"
#include <vnl/vnl_math.h>

namespace itk
{

/**
 * Constructor
 */
PowellOptimizer
::PowellOptimizer()
{
  m_Maximize = false;
  m_ParameterOrderInitialized = false;
  
  m_AcceptNewDirections = true;
  m_NumberOfFunctionEvaluations = 0;
  m_StopCondition = NeverStoppedYet;
  m_MaximumNumberOfIterations = 20;
  m_CurrentIteration = 0;

  m_InitialBracketStepSize = 1.0;
  m_InitialBracketStepSizeShrinkFactor = 1.0;
  m_MaximumBracketStep = 10.0;
  m_FractionalPrecisionBrent = 0.001;
  m_AbsolutePrecisionBrent = 0.0001;
  m_MaximumNumberOfIterationsBrent = 50;
  m_FractionalPrecision = 0.00001;
  m_AbsolutePrecision = 0.0;

}



void
PowellOptimizer
::PrintSelf(std::ostream& os, Indent indent) const
{
  Superclass::PrintSelf(os,indent);

  os << indent << "Maximize: "
       <<  m_Maximize << std::endl;
  os << indent << "ParameterOrderInitialized: "
       <<  m_ParameterOrderInitialized << std::endl;
  os << indent << "AcceptNewDirections: "
       << m_AcceptNewDirections << std::endl;
   os << indent << "ParameterOrder: "
       << m_ParameterOrder << std::endl;
  os << indent << "NumberOfFunctionEvaluations: "
       << m_NumberOfFunctionEvaluations << std::endl;
  os << indent << "CurrentValue: "
      << m_CurrentValue << std::endl;
  os << indent << "StopCondition: "
      << m_StopCondition << std::endl;
  os << indent << "MaximumNumberOfIterations: "
      << m_MaximumNumberOfIterations << std::endl;
  os << indent << "CurrentIteration: "
      << m_CurrentIteration << std::endl;
  os << indent << "CurrentDirection: "
      << m_CurrentDirection << std::endl;
  os << indent << "InitialBracketStepSize: "
      << m_InitialBracketStepSize << std::endl;
  os << indent << "InitialBracketStepSizeShrinkFactor: "
      << m_InitialBracketStepSizeShrinkFactor << std::endl;
  os << indent << "MaximumBracketStep: "
      << m_MaximumBracketStep << std::endl;
  os << indent << "FractionalPrecisionBrent: "
      <<  m_FractionalPrecisionBrent << std::endl;
  os << indent << "AbsolutePrecisionBrent: "
      << m_AbsolutePrecisionBrent << std::endl;
  os << indent << "MaximumNumberOfIterationsBrent: "
      << m_MaximumNumberOfIterationsBrent << std::endl;
  os << indent << "FractionalPrecision: "
      << m_FractionalPrecision << std::endl;
  os << indent << "AbsolutePrecision: "
      << m_AbsolutePrecision << std::endl;

  os << std::endl;

}


/**
 * Set the parameter order
 */
void
PowellOptimizer
::SetParameterOrder(const ParameterOrderType & parameterOrder)
{
  itkDebugMacro("setting parameter order to " <<  parameterOrder);
  m_ParameterOrder = parameterOrder;
  m_ParameterOrderInitialized = true;
  this->Modified();
}


/**
 * Start the optimization
 */
void
PowellOptimizer
::StartOptimization( void )
{

  // Set counter of NumberOfFunctionEvaluations to zero
  m_NumberOfFunctionEvaluations = 0;

  // Constants used by Powell. Here TINY indicates a small number
  const double TINY = 1.0e-25;

  // Set up initial direction table
  const unsigned int spaceDimension =
                        m_CostFunction->GetNumberOfParameters();

  ScalesType scales = this->GetScales();
  std::vector<DirectionType>   originalDirectionTable( spaceDimension );
  for (unsigned int i=0; i<spaceDimension; ++i)
    {
      DirectionType thisDirection( spaceDimension );
      thisDirection.Fill( 0.0 );
      thisDirection[ i ] = 1/scales[ i ];
      originalDirectionTable[ i ] = thisDirection;
    }
    
   // Make sure the parameter order is properly initialized
   if (!m_ParameterOrderInitialized)
    {
      ParameterOrderType parameterOrder( spaceDimension );
      for (unsigned int j=0; j<spaceDimension; j++)
        {
          parameterOrder[ j ] = j+1;
        }
      this->SetParameterOrder( parameterOrder );
    }


  // Check parameter order and decide what parameters need to be optimized, 
  // posssibly bound together with other parameters
  unsigned int   numberOfDirections = 0;
  Array<unsigned int>   lookupTable( spaceDimension );
  lookupTable.Fill( 0 );
  for (unsigned int i=1; i<=spaceDimension; i++)
    {
      bool found = false;
      for (unsigned int j=0; j<spaceDimension; j++)
        {
          if ( m_ParameterOrder[ j ] == i )
            {
              if (!found)
                {
                  numberOfDirections++;
                  found = true;
                }
              lookupTable[ j ] = numberOfDirections;
            }
        }
    }
  //std::cout << "lookupTable: " << lookupTable << std::endl;


  // Set the direction table
  std::vector<DirectionType>   directionTable( numberOfDirections );
  for (unsigned i=0; i<numberOfDirections; i++)
    {
      DirectionType thisDirection( spaceDimension );
      thisDirection.Fill( 0.0 );
      directionTable[ i ] = thisDirection;
    }
  for (unsigned j=0; j<spaceDimension; j++)
    {
      if ( lookupTable[j] )
        {
          for (unsigned k=0; k<spaceDimension; k++)
            {
              directionTable[ lookupTable[j]-1 ] [ k ] += originalDirectionTable[ j ] [ k ];
            }
        }
    }
  //std::cout << "directionTable: " << std::endl;
  for (unsigned int i=0; i<numberOfDirections; ++i)
    {
    //std::cout << "   " <<  directionTable[i] << std::endl;
    }
  //std::cout << std::endl;


  // Initialize some standard stuff
  this->SetCurrentPosition( this->GetInitialPosition() );
  InvokeEvent( StartEvent() );



  // Save the current position
  ParametersType previousPosition = this->GetCurrentPosition();  // pt in NR

  // Get function value at current position
  m_CurrentValue = EvaluateCostFunction( this->GetCurrentPosition() );

  try
    {

      // Start iterations
      for( m_CurrentIteration=0, m_CurrentBracketStepSize = m_InitialBracketStepSize; ; ++m_CurrentIteration,
              m_CurrentBracketStepSize /= m_InitialBracketStepSizeShrinkFactor )
        {
          // Save current function value
          MeasureType previousValue = m_CurrentValue; // fp in NR


          // Optimize over all directions and record which direction gave biggest optimization
          unsigned int biggestFunctionChangeDirectionNumber = 0;  // ibig in NR
          MeasureType biggestFunctionChange = 0.0;  // del in NR

          for (unsigned int i=0; i<numberOfDirections; ++i)
            {
              // Save old value
              MeasureType oldValue = m_CurrentValue;

              // Optimize along current direction
              m_CurrentDirection = directionTable[ i ];
              OptimizeAlongCurrentDirection();

              // Determine if it's the largest change so far
              if ( fabs(oldValue - m_CurrentValue) >  biggestFunctionChange )
                {
                  biggestFunctionChange = vnl_math_abs( oldValue - m_CurrentValue );
                  biggestFunctionChangeDirectionNumber = i;
                }
            }



          // Determine if we should stop, either by convergence or by exceeding maximum
          // number of iterations
          //std::cout << "Testing stop criterion" << std::endl;
          if ( 2.0*vnl_math_abs(previousValue - m_CurrentValue)  <=
                m_FractionalPrecision * ( vnl_math_abs(previousValue) + vnl_math_abs(m_CurrentValue) ) + TINY)
            {
              m_StopCondition = Converged;
              InvokeEvent( EndEvent() );
              //std::cout << "Convergence in function value detected: finished optimizing." << std::endl;
              break;
            }
          if( m_CurrentIteration >= m_MaximumNumberOfIterations-1 )
            {
              m_StopCondition = MaximumNumberOfIterations;
              InvokeEvent( EndEvent() );
              //std::cout << "Maximum number of iterations reached: finished optimizing." << std::endl;
              break;
            }
          ParametersType currentPosition = this->GetCurrentPosition();
          if ( ( element_product( currentPosition - previousPosition, scales ) ).magnitude() <= m_AbsolutePrecision )
            {
              m_StopCondition = Converged;
              InvokeEvent( EndEvent() );
              //std::cout << "Convergence in parameter position detected: finished optimizing." << std::endl;
              break;
            }


          // Construct extrapolated point
          //std::cout << "Constructing extrapolated point" << std::endl;
          ParametersType   extrapolatedPosition( spaceDimension );  // ptt in NR
          DirectionType      extrapolatedDirection( spaceDimension );
          for (unsigned int j=0; j<spaceDimension; j++)
            {
              extrapolatedPosition[j] = 2 * currentPosition[j] - previousPosition[j];
              extrapolatedDirection[j] = currentPosition[j] - previousPosition[j];
            }
          m_CurrentDirection = extrapolatedDirection;
          previousPosition = currentPosition;


          // Evaluate function at extrapolated point
          //std::cout << "Evaluating at extrapolated point" << std::endl;
          MeasureType  extrapolatedValue = SafeEvaluateCostFunction( extrapolatedPosition );

          // Check if the extrapolated direction is promising
          if ( PretendMinimization( extrapolatedValue ) < PretendMinimization( previousValue ) )
            {
              //std::cout << "Extrapolated direction might be promising" << std::endl;
              if ( 2.0*PretendMinimization(  previousValue  - 2.0*m_CurrentValue + extrapolatedValue )
                    * vnl_math_sqr( vnl_math_abs( previousValue - m_CurrentValue ) - biggestFunctionChange )
                     <  biggestFunctionChange * vnl_math_sqr( previousValue - extrapolatedValue ) )
                {
                  //std::cout << "Extrapolated direction is promising" << std::endl;
                  // Move to the optimum of the extrapolated direction
                  OptimizeAlongCurrentDirection();
                  currentPosition = this->GetCurrentPosition();

                  // Save the new direction, and throw away the direction that previously caused the
                  // biggest optimization
                  if ( m_AcceptNewDirections )
                    {
                      //std::cout << "Accepting new direction: " << std::endl;
                      DirectionType newDirection( spaceDimension );
                      for (unsigned int j=0; j<spaceDimension; j++)
                        {
                          newDirection[j] = currentPosition[j] - previousPosition[j];
                        }
                      //std::cout << "     " <<  newDirection << std::endl;
                      //std::cout << "     biggestFunctionChangeDirectionNumber " <<  biggestFunctionChangeDirectionNumber << std::endl;
                      directionTable[ biggestFunctionChangeDirectionNumber ] = directionTable[ numberOfDirections - 1 ];
                      directionTable[ numberOfDirections - 1 ] =  newDirection;

                      //std::cout << "new directionTable: " << std::endl;
                      for (unsigned int i=0; i<numberOfDirections; ++i)
                        {
                        //std::cout << "   " <<  directionTable[i] << std::endl;
                        }
                      //std::cout << std::endl;
                    }
                }
            } // End of inspectation of extrapolated direction


        } // End of loop over iterations

    } // End of try block
  catch( ExceptionObject& err )
    {
      // An exception has occurred.
      InvokeEvent( EndEvent() );
      throw err;
    }

  return;
}



/**
 * Optimize f( m_CurrentPosition + x*m_CurrentDirection ) for scalar x.
 */
void
PowellOptimizer
::OptimizeAlongCurrentDirection( void )
{

  //std::cout << "Optimizing x for "<< this->GetCurrentPosition() << " + x * " << m_CurrentDirection <<  std::endl;


  double ax;
  double bx;
  double cx;

  // Bracket the optimum first
  Bracket(ax, bx, cx);

  // Now close down on optimum
  Brent(ax, bx, cx);

  InvokeEvent( IterationEvent() );

  return;
}



/**
  * Isolate optimum with Brent's method.
  * From Numerical Recipes in C, pp. 404-405, with very few modifications
  */
void
PowellOptimizer
::Brent(double ax, double bx, double cx)
{
  
  //std::cout << "Brent..."<< std::endl;


  /** Constants used by brent. Here CGOLD is the golden ratio; ZEPS is a
  small number that protects against trying to achieve fractional accuracy for
  a minimum that happens to be exactly zero. */
  const double CGOLD = 0.3819660;
  const double ZEPS = 1.0e-10;


  MeasureType fu,fv,fw,fx;
  double a,b,d,etemp, p,q,r,tol1,tol2,u,v,w,x,xm;
  double e=0.0;  // This will be the distance moved on the step before last.

  // a and b must be in ascending order, but input abscissas need not be.
  a=(ax < cx ? ax : cx);
  b=(ax > cx ? ax : cx);

  // Initializations...
  x=w=v=bx;
  fx = EvaluateCostFunction( x );
  fw=fv=fx;
  //std::cout << "   [" << a << "  " << x << "  " << b <<  "]    (" << fx <<  ")"  << std::endl;


  // Main program loop.
  for (unsigned long iter=1; iter<=m_MaximumNumberOfIterationsBrent; iter++)
    {

      // Test for done here.
      xm=0.5*(a+b);
      tol2=2.0*( tol1=m_FractionalPrecisionBrent*vnl_math_abs(x)+ZEPS );
      if ( ( fabs(x-xm) + 0.5*(b-a) ) <= tol2  ||
           ( fabs(x-xm) + 0.5*(b-a) ) <=  m_AbsolutePrecisionBrent )
        {
          const unsigned int spaceDimension =
                        m_CostFunction->GetNumberOfParameters();
          ParametersType currentPosition = this->GetCurrentPosition();
          ParametersType newPosition( spaceDimension );
          for (unsigned int j=0; j<spaceDimension; j++)
            {
              newPosition[j] = currentPosition[j] + x*m_CurrentDirection[j];
            }
          this->SetCurrentPosition( newPosition );
          m_CurrentValue =  fx;
          return;
        }
      if (fabs(e) > tol1)
        {
          // Construct a trial parabolic fit.
          r=(x-w)*PretendMinimization(fx-fv);
          q=(x-v)*PretendMinimization(fx-fw);
          p=(x-v)*q-(x-w)*r;
          q=2.0*(q-r);
          if (q > 0.0) p = -p;
          q=fabs(q);
          etemp=e;
          e=d;
          if (fabs(p) >= fabs(0.5*q*etemp) || p <= q*(a-x) || p >= q*(b-x))
            {
              d=CGOLD*(e=(x >= xm ? a-x : b-x));
            }
          // The above conditions determine the acceptability of the parabolic fit. Here
          // we take the golden section step into the larger of the two segments.
          else
            {
              // Take the parabolic step.
              d=p/q;
              u=x+d;
              if (u-a < tol2 || b-u < tol2)
                {
                  d = ( (xm-x) > 0 ? vnl_math_abs(tol1) : -vnl_math_abs(tol1) );
                }
            }
        }
      else
        {
          d=CGOLD*(e=(x >= xm ? a-x : b-x));
        }

      double tmp = ( (d)>0 ? vnl_math_abs(tol1) : -vnl_math_abs(tol1) );
      u=(vnl_math_abs(d) >= tol1 ? x+d : x+tmp);
      fu = EvaluateCostFunction( u );      // This is the one function evaluation per iteration.

      // Now decide what to do with our function evaluation.
      // Housekeeping follows:
      if (PretendMinimization(fu) <= PretendMinimization(fx))
        {
          if (u >= x) a=x; else b=x;
          v = w;
          w = x;
          x = u;
          fv = fw;
          fw = fx;
          fx = fu;
          //std::cout << "   [" << a << "  " << x << "  " << b <<  "]    (" << fx <<  ")"  << std::endl;
        }
      else
        {
          if (u < x) a=u; else b=u;
          //std::cout << "   [" << a << "  " << x << "  " << b <<  "]    (" << fx <<  ")"  << std::endl;
          
          if (PretendMinimization(fu) <= PretendMinimization(fw) || w == x)
            {
              v=w;
              w=u;
              fv=fw;
              fw=fu;
            }
          else if (PretendMinimization(fu) <= PretendMinimization(fv) || v == x || v == w)
            {
              v=u;
              fv=fu;
            }
        }

    } // Done with housekeeping. Back for another iteration.

  // Reached maximum number of iterations
  // Go to where we got so far
  const unsigned int spaceDimension =
        m_CostFunction->GetNumberOfParameters();
  ParametersType currentPosition = this->GetCurrentPosition();
  ParametersType newPosition( spaceDimension );
  for (unsigned int j=0; j<spaceDimension; j++)
    {
      newPosition[j] = currentPosition[j] + x*m_CurrentDirection[j];
    }
  this->SetCurrentPosition( newPosition );
  m_CurrentValue =  fx;

  return;
}



/**
 * Bracket the optimum
 * From Numerical Recipes in C, pp. 400-401, with very few modifications
 */
void
PowellOptimizer
::Bracket( double& ax, double& bx, double& cx )
{

  //std::cout << "Bracketing..."<< std::endl;

  /** GOLD is the default ratio by which successive
  intervals are magnified. TINY is used to prevent any
  possible division by zero.*/
  const double   GOLD = 1.618034;
  const double   TINY = 1.0e-20;


  MeasureType fa, fb, fc, fdum;
  double ulim,u,r,q,fu,dum;

  // Initialization
  ax = 0.0;
  fa = m_CurrentValue;
  bx = m_CurrentBracketStepSize;
  fb = SafeEvaluateCostFunction( bx );

  if (PretendMinimization(fb) > PretendMinimization(fa))
    {
      // Switch roles of a and b so that we can go downhill in the direction from a to b.
      dum = ax;
      ax = bx;
      bx = dum;
      fdum = fb;
      fb = fa;
      fa = fdum;
    }

  cx=(bx)+GOLD*(bx-ax); // First guess for c.
  fc = SafeEvaluateCostFunction( cx );
  //std::cout << "   [" << ax << "  " << bx << "  " << cx <<  "]    (" << fa << "  " << fb << "  " << fc << ")"  << std::endl;

  while (PretendMinimization(fb) > PretendMinimization(fc))
    {
      // Keep returning here until we bracket.

      // Compute u by parabolic extrapolation from a, b, c. TINY is used to prevent any
      // possible division by zero.
      r=(bx-ax)*PretendMinimization(fb-fc);
      q=(bx-cx)*PretendMinimization(fb-fa);
      double tmp = ( vnl_math_abs(q-r)>TINY ? vnl_math_abs(q-r) : TINY );
      double tmp2 = ( (q-r)>0.0 ? vnl_math_abs(tmp) : -vnl_math_abs(tmp) );
      u=bx - ((bx-cx)*q-(bx-ax)*r)/(2.0*tmp2);
      ulim=bx + m_MaximumBracketStep*(cx-bx); // We won't go farther than this.

      // Test various possibilities:
      if ((bx-u)*(u-cx) > 0.0)
        {
          // Parabolic u is between b and c: try it.
          fu = SafeEvaluateCostFunction( u );
          if (PretendMinimization(fu) < PretendMinimization(fc))
            {
              // Got a minimum between b and c.
              ax=bx;
              bx=u;
              fa=fb;
              fb=fu;
              //std::cout << "   [" << ax << "  " << bx << "  " << cx <<  "]    (" << fa << "  " << fb << "  " << fc << ")"  << std::endl;
              return;
            }
          else if (PretendMinimization(fu) > PretendMinimization(fb))
            {
              // Got a minimum between between a and u.
              cx=u;
              fc=fu;
              //std::cout << "   [" << ax << "  " << bx << "  " << cx <<  "]    (" << fa << "  " << fb << "  " << fc << ")"  << std::endl;
              return;
            }

          //  Parabolic fit was no use. Use default magnification.
          u=cx+GOLD*(cx-bx);
          fu = SafeEvaluateCostFunction( u );
        }
      else if ((cx-u)*(u-ulim) > 0.0)
        {
          // Parabolic fit is between c and its allowed limit.
          fu = SafeEvaluateCostFunction( u );
          if (PretendMinimization(fu) < PretendMinimization(fc))
            {
              bx = cx;
              cx = u;
              u = cx+GOLD*(cx-bx);
              fb = fc;
              fc = fu;
              fu = SafeEvaluateCostFunction( u );
              //std::cout << "   [" << ax << "  " << bx << "  " << cx <<  "]    (" << fa << "  " << fb << "  " << fc << ")"  << std::endl;
           }
        }
      else if ((u-ulim)*(ulim-cx) >= 0.0)
        {
          // Limit parabolic u to maximum allowed value.
          u=ulim;
          fu = SafeEvaluateCostFunction( u );
        }
      else
        {
          // Reject parabolic u, use default magnification.
          u=cx+GOLD*(cx-bx);
          fu = SafeEvaluateCostFunction( u );
        }

      // Eliminate oldest point and continue.
      ax = bx;
      bx = cx;
      cx = u;
      fa = fb;
      fb = fc;
      fc = fu;
      //std::cout << "   [" << ax << "  " << bx << "  " << cx <<  "]    (" << fa << "  " << fb << "  " << fc << ")"  << std::endl;
    }

}




} // end namespace itk

#endif
