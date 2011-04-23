#ifndef __itkPowellOptimizer_h
#define __itkPowellOptimizer_h

#include "itkSingleValuedNonLinearOptimizer.h"

namespace itk
{


class PowellOptimizer : public SingleValuedNonLinearOptimizer
{
public:
  /** Standard class typedefs. */
  typedef PowellOptimizer          Self;
  typedef SingleValuedNonLinearOptimizer    Superclass;
  typedef SmartPointer<Self>                Pointer;
  typedef SmartPointer<const Self>          ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Run-time type information (and related methods). */
  itkTypeMacro( PowellOptimizer, SingleValuedNonLinearOptimizer );

  /**  Parameters type.
   *  It defines a position in the optimization search space. */
  typedef Superclass::ParametersType   ParametersType;

  /** */
  typedef Superclass::ParametersType   DirectionType;

  /** */
  typedef Array<unsigned int>              ParameterOrderType;

  /**  Measure type.
   *  It defines a type used to return the cost function value.  */
  typedef Superclass::MeasureType   MeasureType;

  /** Codes of stopping conditions */
  typedef enum {
    NeverStoppedYet,
    Converged,
    MaximumNumberOfIterations,
    MetricError
  } StopConditionType;

  /** Methods to configure the cost function. */
  itkGetMacro( Maximize, bool );
  itkSetMacro( Maximize, bool );
  itkBooleanMacro( Maximize );
  bool GetMinimize( ) const
    { return !m_Maximize; }
  void SetMinimize(bool v)
    { this->SetMaximize(!v); }
  void MinimizeOn()
    { this->MaximizeOff(); }
  void MinimizeOff()
    { this->MaximizeOn(); }

  /** Start optimization. */
  void    StartOptimization( void );

  /** Set the number of iterations. */
  itkSetMacro( MaximumNumberOfIterations, unsigned long );

  /** Get the number of iterations. */
  itkGetConstMacro( MaximumNumberOfIterations, unsigned long );

  /** Get the current iteration number. */
  itkGetConstMacro( CurrentIteration, unsigned int );

  /** Get the current value. */
  itkGetConstMacro( CurrentValue, double );

  /** Get Stop condition. */
  itkGetConstMacro( StopCondition, StopConditionType );

  /** */
  itkGetMacro( AcceptNewDirections, bool );
  itkSetMacro( AcceptNewDirections, bool );
  itkBooleanMacro( AcceptNewDirections );

  /** */
  itkGetConstMacro( NumberOfFunctionEvaluations, unsigned long );

  /** */
  itkGetConstMacro( CurrentDirection, DirectionType );

  /** */
  itkSetMacro( InitialBracketStepSize,  double );

  /** */
  itkGetConstMacro( InitialBracketStepSize,  double );

  /**  */
  itkSetMacro( InitialBracketStepSizeShrinkFactor, double );

  /** */
  itkGetConstMacro( InitialBracketStepSizeShrinkFactor,  double );

  /**  */
  itkSetMacro( MaximumBracketStep, double );

  /** */
  itkGetConstMacro( MaximumBracketStep, double );

  /**  */
  itkSetMacro( FractionalPrecisionBrent, double );

  /** */
  itkGetConstMacro( FractionalPrecisionBrent, double );

  /**  */
  itkSetMacro( AbsolutePrecisionBrent, double );

  /** */
  itkGetConstMacro( AbsolutePrecisionBrent, double );

  /**  */
  itkSetMacro( MaximumNumberOfIterationsBrent,  unsigned long );

  /** */
  itkGetConstMacro( MaximumNumberOfIterationsBrent, unsigned long );

  /**  */
  itkSetMacro( FractionalPrecision,  double );

  /** */
  itkGetConstMacro( FractionalPrecision,  double );

  /**  */
  itkSetMacro( AbsolutePrecision,  double );

  /** */
  itkGetConstMacro( AbsolutePrecision,  double );


  /** Set parameters order. */
  void SetParameterOrder( const ParameterOrderType & parameterOrder );

  /** Get current parameters order. */
  itkGetConstMacro( ParameterOrder, ParameterOrderType );




protected:
  PowellOptimizer();
  virtual ~PowellOptimizer() {};
  void PrintSelf(std::ostream& os, Indent indent) const;


  // made protected so subclass can access
  bool                           m_Maximize;
  bool                           m_ParameterOrderInitialized;

private:
  PowellOptimizer(const Self&); //purposely not implemented
  void operator=(const Self&); //purposely not implemented


  /**
  Given distinct initial points ax and bx, this routine
  searches in the downhill direction (defined by the function as evaluated at the
  initial points) and returns new points ax, bx, cx that bracket a minimum of the
  function.
  */
  void Bracket(double& ax, double& bx, double& cx);


  /**
  Given a bracketing triplet of abscissas ax, bx, cx
  (such that bx is between ax and cx,and f(bx) is better than both f(ax) and f(cx)),
  this routine isolates the optimum using Brent's method. Subsequently,
  m_CurrentPosition is updated accordingly, and m_CurrentValue is set to the function
  value corresponding to the new position.
  */
  void Brent(double ax, double bx, double cx);


  /**
  Optimizes m_CurrentPosition + x*m_CurrentDirection for scalar x. Subsequently,
  m_CurrentPosition is updated accordingly, and m_CurrentValue is set to the function
  value corresponding to the new position.
  */
  void  OptimizeAlongCurrentDirection( void );


  /** If we are maximizing, switch the sign of the function so that Powell and the line mimization
  "think" we are actually minimizing. */
  inline MeasureType PretendMinimization(const MeasureType& value) const
    { return ( m_Maximize ? -value : value ); }


  /**
  Evaluate cost function at m_CurrentPosition + x*m_CurrentDirection
  */
  MeasureType EvaluateCostFunction( double x )
    {
    const ParametersType & currentPosition = this->GetCurrentPosition();

    const unsigned int spaceDimension =
                        m_CostFunction->GetNumberOfParameters();
    ParametersType position( spaceDimension );
    for (unsigned int j=0; j<spaceDimension; j++)
      {
        position[j] = currentPosition[j] + x * m_CurrentDirection[j];
      }

    return EvaluateCostFunction( position );
    }


  /**
  Evaluate cost function at m_CurrentPosition + x*m_CurrentDirection
  */
  MeasureType SafeEvaluateCostFunction( double x )
    {
    const ParametersType & currentPosition = this->GetCurrentPosition();

    const unsigned int spaceDimension =
                        m_CostFunction->GetNumberOfParameters();
    ParametersType position( spaceDimension );
    for (unsigned int j=0; j<spaceDimension; j++)
      {
        position[j] = currentPosition[j] + x * m_CurrentDirection[j];
      }

    return SafeEvaluateCostFunction( position );
    }



  /**
  Evaluate cost function at a specified position
  */
  MeasureType EvaluateCostFunction( const ParametersType& position )
    {

      MeasureType value;
      try
        {
        value = m_CostFunction->GetValue( position );
        }
      catch( ExceptionObject& err )
        {
          // An exception has occurred.
          m_StopCondition = MetricError;
          throw err;
        }

      m_NumberOfFunctionEvaluations++;

      return value;
    }


  // Evaluates cost function at specified position, but returns a very
  // bad value if an exception occurs during the evaluation.
  MeasureType SafeEvaluateCostFunction( const ParametersType& position )
    {
    MeasureType value;

    try
      {
      value = m_CostFunction->GetValue( position );
      }
    catch( ExceptionObject& )
      {
      // An exception has occurred.
      std::cout << "Got exception from the metric; converting to very bad value" << std::endl;
      if ( m_Maximize )
        {
        // Lowest possible value
        value = itk::NumericTraits<MeasureType>::min();
        }
      else
        {
        // Highest possible value
        value = itk::NumericTraits<MeasureType>::max();
        }
     }

    m_NumberOfFunctionEvaluations++;

    return value;
    }


  bool                               m_AcceptNewDirections;
  ParameterOrderType        m_ParameterOrder;
  unsigned long                 m_NumberOfFunctionEvaluations;
  MeasureType                 m_CurrentValue;
  StopConditionType          m_StopCondition;
  unsigned long                 m_MaximumNumberOfIterations;
  unsigned long                 m_CurrentIteration;
  DirectionType                 m_CurrentDirection;
  double                           m_InitialBracketStepSize;  // the value 1.0 in routine "linmin" in
                                                                           // Numerical Recipes in C, pp. 419
  double                           m_InitialBracketStepSizeShrinkFactor; // How m_InitialBracketStepSize
                                                                                            // is decreased over Powell's iterations
  double                           m_MaximumBracketStep;  // GLIMIT in routine "mnbrak" in
                                                                           // Numerical Recipes in C, pp. 400-401
  double                           m_FractionalPrecisionBrent;  // tol in routine "brent" in
                                                                               // Numerical Recipes in C, pp. 404-405
  double                           m_AbsolutePrecisionBrent;  // Absolute precision with which x is isolated
                                                                              // Here m_CurrentPosition + x*m_CurrentDirection
  unsigned long                 m_MaximumNumberOfIterationsBrent; // ITMAX in routine "brent" in
                                                                                           // Numerical Recipes in C, pp. 404-405
  double                           m_FractionalPrecision;  // Fractional tolerance in cost function value used for
                                                                        // determing convergence of Powell's algorithm. This is 
                                                                        // called ftol in routine "powell" in
                                                                        // Numerical Recipes in C, pp. 417-418
  double                           m_AbsolutePrecision;

  double                           m_CurrentBracketStepSize;


};

} // end namespace itk


#endif



