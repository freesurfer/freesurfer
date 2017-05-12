#include "kvlMatlabRunner.h" 
#include "kvlMatlabObjectArray.h"
#include "kvlAtlasMesh.h"
#include "kvlAtlasMeshDeformationOptimizer.h"


namespace kvl
{

  
class StepOptimizer : public MatlabRunner
{
public:
  /** Smart pointer typedef support. */
  typedef StepOptimizer         Self;
  typedef itk::Object              Superclass;
  typedef itk::SmartPointer<Self>  Pointer;
  typedef itk::SmartPointer<const Self>  ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro( Self );

  /** Run-time type information (and related methods). */
  itkTypeMacro( StepOptimizer, itk::Object );

  virtual void Run( int nlhs, mxArray* plhs[],
                    int nrhs, const mxArray* prhs[] )
    {
    //std::cout << "I am " << this->GetNameOfClass() 
    //          << " and I'm running! " << std::endl;
              
              
    // [ minLogLikelihoodTimesPrior, maximalDeformation ] = kvlStepOptimizer( optimizer )
  
    // Make sure input arguments are correct
    if ( ( nrhs < 1 ) || 
         !mxIsInt64( prhs[ 0 ] ) )
      {
      mexErrMsgTxt( "Incorrect arguments" );
      }

    // Retrieve the optimizer and let it do the work
    const int optimizerHandle = *( static_cast< int* >( mxGetData( prhs[ 0 ] ) ) );
    itk::Object::ConstPointer object = kvl::MatlabObjectArray::GetInstance()->GetObject( optimizerHandle );
    kvl::AtlasMeshDeformationOptimizer::ConstPointer  constOptimizer
             = dynamic_cast< const kvl::AtlasMeshDeformationOptimizer* >( object.GetPointer() );
    if ( !constOptimizer.GetPointer() )
      {
      std::cout << "typeid: " << typeid( *object ).name() << std::endl;
      mexErrMsgTxt( "optimizer doesn't refer to the correct ITK object type" );
      }  
    kvl::AtlasMeshDeformationOptimizer::Pointer  optimizer
                   = const_cast< kvl::AtlasMeshDeformationOptimizer* >( constOptimizer.GetPointer() );
    

    // Let the beast go
    const double  maximalDeformation = optimizer->Step();
    const double  minLogLikelihoodTimesPrior = optimizer->GetMinLogLikelihoodTimesPrior();
    
    
    // Return the maximalDeformation to Matlab
    plhs[ 0 ] = mxCreateDoubleScalar( minLogLikelihoodTimesPrior );
    plhs[ 1 ] = mxCreateDoubleScalar( maximalDeformation );
    }
  
protected:
  StepOptimizer() {};
  virtual ~StepOptimizer() {};


  StepOptimizer(const Self&); //purposely not implemented
  void operator=(const Self&); //purposely not implemented

private:

};

} // end namespace kvl



