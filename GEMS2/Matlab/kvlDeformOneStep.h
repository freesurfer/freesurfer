#include "kvlMatlabRunner.h" 
#include "kvlMatlabObjectArray.h"
#include "kvlAtlasMesh.h"
#include "kvlCroppedImageReader.h"
#include "kvlAtlasMeshDeformationLevenbergMarquardtOptimizer.h"
#include "kvlAtlasMeshDeformationConjugateGradientOptimizer.h"


namespace kvl
{

  
class DeformOneStep : public MatlabRunner
{
public:
  /** Smart pointer typedef support. */
  typedef DeformOneStep         Self;
  typedef itk::Object              Superclass;
  typedef itk::SmartPointer<Self>  Pointer;
  typedef itk::SmartPointer<const Self>  ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro( Self );

  /** Run-time type information (and related methods). */
  itkTypeMacro( DeformOneStep, itk::Object );

  virtual void Run( int nlhs, mxArray* plhs[],
                    int nrhs, const mxArray* prhs[] )
    {
    //std::cout << "I am " << this->GetNameOfClass() 
    //          << " and I'm running! " << std::endl;
              
              
    // [ minLogLikelihoodTimesPrior, maximalDeformation ] = kvlDeformOneStep( optimizer )
  
    // Make sure input arguments are correct
    if ( ( nrhs < 1 ) || 
         !mxIsInt64( prhs[ 0 ] ) )
      {
      mexErrMsgTxt( "Incorrect arguments" );
      }

    // Retrieve the optimizer and let it do the work
    double  maximalDeformation = 0.0;
    double  minLogLikelihoodTimesPrior = 0.0;
    const int optimizerHandle = *( static_cast< int* >( mxGetData( prhs[ 0 ] ) ) );
    itk::Object::ConstPointer object = kvl::MatlabObjectArray::GetInstance()->GetObject( optimizerHandle );
    if ( typeid( *object ) == typeid( AtlasMeshDeformationLevenbergMarquardtOptimizer ) )
      {
      AtlasMeshDeformationLevenbergMarquardtOptimizer::ConstPointer constOptimizer 
          = static_cast< const AtlasMeshDeformationLevenbergMarquardtOptimizer* >( object.GetPointer() );
      AtlasMeshDeformationLevenbergMarquardtOptimizer::Pointer  optimizer 
          = const_cast< AtlasMeshDeformationLevenbergMarquardtOptimizer* >( constOptimizer.GetPointer() );

      // Let the beast go
      maximalDeformation = optimizer->PerformOneSuccessfulStep();
      minLogLikelihoodTimesPrior = optimizer->GetMinLogLikelihoodTimesPrior();
      }
    else if ( typeid( *object ) == typeid( AtlasMeshDeformationConjugateGradientOptimizer ) )
      {
      AtlasMeshDeformationConjugateGradientOptimizer::ConstPointer constOptimizer 
          = static_cast< const AtlasMeshDeformationConjugateGradientOptimizer* >( object.GetPointer() );
      AtlasMeshDeformationConjugateGradientOptimizer::Pointer  optimizer 
          = const_cast< AtlasMeshDeformationConjugateGradientOptimizer* >( constOptimizer.GetPointer() );

      // Let the beast go
      maximalDeformation = optimizer->PerformOneIteration();
      minLogLikelihoodTimesPrior = optimizer->GetMinLogLikelihoodTimesPrior();
      }
    else
      {  
      mexErrMsgTxt( "mesh doesn't refer to the correct ITK object type" );
      }
    
    
    // Return the maximalDeformation to Matlab
    plhs[ 0 ] = mxCreateDoubleScalar( minLogLikelihoodTimesPrior );
    plhs[ 1 ] = mxCreateDoubleScalar( maximalDeformation );
    }
  
protected:
  DeformOneStep() {};
  virtual ~DeformOneStep() {};


  DeformOneStep(const Self&); //purposely not implemented
  void operator=(const Self&); //purposely not implemented

private:

};

} // end namespace kvl



