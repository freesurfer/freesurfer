#ifndef __kvlParameterOrderPowellOptimizer_h
#define __kvlParameterOrderPowellOptimizer_h

#include "itkPowellOptimizer.h"
#include "itkSingleValuedCostFunction.h"
#include "itkCommand.h"


namespace kvl
{


class WrappedSingleValuedCostFunction : public itk::SingleValuedCostFunction
{

public:
  /** Standard class typedefs. */
  typedef WrappedSingleValuedCostFunction     Self;
  typedef itk::SingleValuedCostFunction                 Superclass;
  typedef itk::SmartPointer<Self>           Pointer;
  typedef itk::SmartPointer<const Self>     ConstPointer;

  /** Run-time type information (and related methods). */
  itkTypeMacro( WrappedSingleValuedCostFunction, itk::SingleValuedCostFunction );

  /** */
  itkNewMacro(Self);

  // Some typedefs
  typedef Superclass::ParametersType  ParametersType;
  typedef Superclass::MeasureType     MeasureType;
  typedef Superclass::DerivativeType  DerivativeType;
  typedef itk::Array< unsigned int >  ParameterOrderType;


  //
  virtual MeasureType GetValue( const ParametersType & parameters ) const
  {
    // Translate the (short) parameters into (long) parameters the cost function understands
    ParametersType    costFunctionParameters = m_DefaultCostFunctionParameters;
    for ( unsigned int i = 0; i < m_ParameterOrder.Size(); i++ )
    {
      if ( m_ParameterOrder[ i ] )
      {
        costFunctionParameters[ i ] = parameters[ m_ParameterOrder[ i ] - 1 ];
      }
    }

    return m_CostFunction->GetValue( costFunctionParameters );
  }

  //
  virtual void GetDerivative( const ParametersType & parameters,
                              DerivativeType & derivative ) const
  {
  }

  virtual unsigned int GetNumberOfParameters() const
  {
    unsigned int  numberOfParameters = 0;
    for ( unsigned int i = 0; i < m_ParameterOrder.Size(); i++ )
    {
      if ( m_ParameterOrder[ i ] != 0 )
      {
        numberOfParameters++;
      }
    }

    return numberOfParameters;
  }





  //
  void  SetUp( itk::SingleValuedCostFunction* costFunction,
               const ParametersType& defaultCostFunctionParameters,
               const ParameterOrderType& parameterOrder )
  {
    m_CostFunction = costFunction;
    m_DefaultCostFunctionParameters = defaultCostFunctionParameters;
    m_ParameterOrder = parameterOrder;
  }




protected:
  WrappedSingleValuedCostFunction()
  {
    m_CostFunction = 0;
  }

  virtual ~WrappedSingleValuedCostFunction()
  {
  }

private:
  WrappedSingleValuedCostFunction(const Self&); //purposely not implemented
  void operator=(const Self&); //purposely not implemented

  itk::SingleValuedCostFunction::Pointer  m_CostFunction;
  ParametersType  m_DefaultCostFunctionParameters;
  ParameterOrderType  m_ParameterOrder;
};





class ParameterOrderPowellOptimizer : public itk::SingleValuedNonLinearOptimizer
{
public:
  /** Standard class typedefs. */
  typedef ParameterOrderPowellOptimizer          Self;
  typedef itk::SingleValuedNonLinearOptimizer    Superclass;
  typedef itk::SmartPointer<Self>                Pointer;
  typedef itk::SmartPointer<const Self>          ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Run-time type information (and related methods). */
  itkTypeMacro( ParameterOrderPowellOptimizer, itk::SingleValuedNonLinearOptimizer );

  // Some typedefs
  typedef Superclass::MeasureType   MeasureType;
  typedef WrappedSingleValuedCostFunction::ParameterOrderType  ParameterOrderType;
  typedef Superclass::ParametersType ParametersType;
  typedef Superclass::ScalesType     ScalesType;

  //
  const ParametersType&  GetCurrentPosition() const
  {
    // Make sure we have what it takes to compute the current position
    if ( ( m_ParameterOrder.Size() == 0 ) ||
         ( m_WrappedOptimizer->GetCurrentPosition().Size() == 0 ) ||
         ( this->GetInitialPosition().Size() == 0 ) )
    {
      return Superclass::GetCurrentPosition();
    }

    // Translate the (short) parameters into (long) parameters the cost function understands
    ParametersType  currentPosition = this->GetInitialPosition();
    for ( unsigned int i = 0; i < m_ParameterOrder.Size(); i++ )
    {
      if ( m_ParameterOrder[ i ] )
      {
        currentPosition[ i ] = m_WrappedOptimizer->GetCurrentPosition()[ m_ParameterOrder[ i ] - 1 ];
      }
    }

    const_cast< Self* >( this )->SetCurrentPosition( currentPosition );
    return Superclass::GetCurrentPosition();
  }

  const std::string GetStopConditionDescription() const
  {
    return m_WrappedOptimizer->GetStopConditionDescription();
  }

  // Start optimization
  void StartOptimization( void )
  {
    // Set up cost function
    WrappedSingleValuedCostFunction::Pointer  metric = WrappedSingleValuedCostFunction::New();
    metric->SetUp( const_cast< CostFunctionType* >( this->GetCostFunction() ), this->GetInitialPosition(), m_ParameterOrder );

    // Set up initial position and scales
    ParametersType  initialPosition( metric->GetNumberOfParameters() );
    ScalesType  scales( metric->GetNumberOfParameters() );
    for ( unsigned int i = 0; i < m_ParameterOrder.Size(); i++ )
    {
      if ( m_ParameterOrder[ i ] )
      {
        scales[ m_ParameterOrder[ i ] - 1 ] = this->GetScales()[ i ];
        initialPosition[ m_ParameterOrder[ i ] - 1 ] = this->GetInitialPosition()[ i ];
      }
    }

    // Set up real Powell Optimizer that works in lower dimensional space
    m_WrappedOptimizer->SetCostFunction( metric );
    m_WrappedOptimizer->SetInitialPosition( initialPosition );
    m_WrappedOptimizer->SetScales( scales );

    // TODO: forward itk::Events from real optimizer to your wrapper optimizer

    // Start optimizer
    m_WrappedOptimizer->StartOptimization();
  }

  // Set parameters order
  void SetParameterOrder( const ParameterOrderType & parameterOrder )
  {
    m_ParameterOrder = parameterOrder;
  }

  const ParameterOrderType& GetParameterOrder() const
  {
    return m_ParameterOrder;
  }



  // Pass on optimizer stuff to the underlying wrapped optimizer
  virtual const unsigned int&  GetCurrentIteration() const
  {
    return m_WrappedOptimizer->GetCurrentIteration();
  }
  virtual void  SetMaximize(bool maximize )
  {
    m_WrappedOptimizer->SetMaximize( maximize );
  }
  virtual void  MaximizeOn()
  {
    m_WrappedOptimizer->SetMaximize( true );
  }
  virtual void  MaximizeOff()
  {
    m_WrappedOptimizer->SetMaximize( false );
  }
  virtual const bool&  GetMaximize() const
  {
    return m_WrappedOptimizer->GetMaximize();
  }
  virtual void  SetMaximumIteration( unsigned int maximumIteration )
  {
    m_WrappedOptimizer->SetMaximumIteration( maximumIteration );
  }
  virtual const unsigned int&  GetMaximumIteration() const
  {
    return m_WrappedOptimizer->GetMaximumIteration();
  }
  virtual void  SetMaximumLineIteration( unsigned int maximumLineIteration )
  {
    m_WrappedOptimizer->SetMaximumLineIteration( maximumLineIteration );
  }
  virtual unsigned int  GetMaximumLineIteration() const
  {
    return m_WrappedOptimizer->GetMaximumLineIteration();
  }
  virtual void  SetStepLength( double stepLength )
  {
    m_WrappedOptimizer->SetStepLength( stepLength );
  }
  virtual const double&  GetStepLength() const
  {
    return m_WrappedOptimizer->GetStepLength();
  }
  virtual void  SetStepTolerance( double stepTolerance )
  {
    m_WrappedOptimizer->SetStepTolerance( stepTolerance );
  }
  virtual const double&  GetStepTolerance() const
  {
    return m_WrappedOptimizer->GetStepTolerance();
  }
  virtual void  SetValueTolerance( double valueTolerance )
  {
    m_WrappedOptimizer->SetValueTolerance( valueTolerance );
  }
  virtual const double&  GetValueTolerance() const
  {
    return m_WrappedOptimizer->GetValueTolerance();
  }
  virtual const MeasureType&   GetCurrentCost() const
  {
    return m_WrappedOptimizer->GetCurrentCost();
  }



  // Forward events from wrapped optimizer to self
  void HandleWrappedOptimizerEvent( itk::Object* object, const itk::EventObject& event )
  {
    this->InvokeEvent( event );
  }


protected:
  ParameterOrderPowellOptimizer()
  {
    //
    m_WrappedOptimizer = itk::PowellOptimizer::New();


    // Add observers to the wrapped optimizer
    typedef itk::MemberCommand< Self >   MemberCommandType;
    MemberCommandType::Pointer  command = MemberCommandType::New();
    command->SetCallbackFunction( this, &Self::HandleWrappedOptimizerEvent );
    m_WrappedOptimizer->AddObserver( itk::StartEvent(), command );
    m_WrappedOptimizer->AddObserver( itk::IterationEvent(), command );
    m_WrappedOptimizer->AddObserver( itk::EndEvent(), command );
    
    // Set correct options to replicate Eugenio's desires
    m_WrappedOptimizer->SetCatchGetValueException( true );
    m_WrappedOptimizer->SetMetricWorstPossibleValue( itk::NumericTraits< double >::max() );

  }

  virtual ~ParameterOrderPowellOptimizer() {};

private:
  ParameterOrderPowellOptimizer(const Self&); //purposely not implemented
  void operator=(const Self&); //purposely not implemented

  ParameterOrderType  m_ParameterOrder;
  itk::PowellOptimizer::Pointer  m_WrappedOptimizer;


};

} // end namespace kvl


#endif



