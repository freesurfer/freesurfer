#include "kvlMatlabRunner.h" 
#include "kvlMatlabObjectArray.h"
#include "kvlAtlasMesh.h"
#include "kvlAtlasMeshDeformationConjugateGradientOptimizerCPU.h"
#include "kvlCroppedImageReader.h"


namespace kvl
{

  
class SetOptimizerPropertiesCPU : public MatlabRunner
{
public:
  /** Smart pointer typedef support. */
  typedef SetOptimizerPropertiesCPU         Self;
  typedef itk::Object              Superclass;
  typedef itk::SmartPointer<Self>  Pointer;
  typedef itk::SmartPointer<const Self>  ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro( Self );

  /** Run-time type information (and related methods). */
  itkTypeMacro( SetOptimizerPropertiesCPU, itk::Object );

  virtual void Run( int nlhs, mxArray* plhs[],
                    int nrhs, const mxArray* prhs[] )
    {
    //std::cout << "I am " << this->GetNameOfClass() 
    //          << " and I'm running! " << std::endl;
              
              
    // kvlSetOptimizerProperties( optimizer, means, variances )
  
    // Make sure input arguments are correct
    if ( ( nrhs < 3 ) || 
         !mxIsInt64( prhs[ 0 ] ) ||
         !mxIsDouble( prhs[ 1 ] ) ||
         !mxIsDouble( prhs[ 2 ] ) )
      {
      mexErrMsgTxt( "Incorrect arguments" );
      }
      
        
    // Retrieve means and variances
    const int  numberOfClasses = mxGetDimensions( prhs[ 1 ] )[ 0 ];
    const int  numberOfImages  = mxGetDimensions( prhs[ 1 ] )[ 1 ];
    mexPrintf("numberOfClasses = %d\n",numberOfClasses);
    mexPrintf("numberOfImages = %d\n",numberOfImages);
    std::vector< vnl_vector< float > > means;
    std::vector< vnl_matrix< float > > precisions;
    vnl_vector< float >  mean ( numberOfImages, 0.0f );
    vnl_matrix< float >  precision( numberOfImages, numberOfImages, 0.0f);

    for ( int classNumber = 0; classNumber < numberOfClasses; classNumber++ )
      {
	for ( int nima = 0; nima < numberOfImages; nima++ )
	  {
	    mean[ nima ] = (mxGetPr( prhs[ 1 ] ))[ classNumber + numberOfClasses*nima ];
	   
	  }
	means.push_back(mean);
      }
    
    //Does not really matter which way you read these in, because the precisions are symmetric matrices
    //transpose wont do any harm. 
    for ( unsigned int classNumber = 0; classNumber < numberOfClasses; classNumber++ )
      {
    	for ( unsigned int row = 0; row < numberOfImages; row++ )
    	  {
    	    for ( unsigned int col = 0; col < numberOfImages; col++ )
    	      {
    		precision[ row ][ col ] = mxGetPr( prhs[ 2 ] )[ row + numberOfImages*(col + numberOfImages*classNumber) ];
    	      }
   	  }
    	precisions.push_back(precision);
      }

    // Show what we have so far

    for ( unsigned int classNumber = 0; classNumber < numberOfClasses; classNumber++ )
      {
	vnl_matrix<float> precMat = precisions[classNumber];
	for ( unsigned int row = 0; row < numberOfImages; row++ )
	  {
	    for ( unsigned int col = 0; col < numberOfImages; col++ )
	      {
		mexPrintf("precisions[%d][%d][%d] = %f\n",row,col,classNumber,precMat[row][col]);
	      }
	  }
      }

    for (unsigned int classNumber = 0; classNumber < numberOfClasses; classNumber++ )
      {
	vnl_vector<float> miini = means[classNumber];
	for ( unsigned int nima = 0; nima < numberOfImages; nima++ )
	  {
	    mexPrintf("means[%d][%d] = %f\n",nima, classNumber, miini[nima]);
	  }
       }
       
    
    // Retrieve the optimizer, and set the means and variances
    const int optimizerHandle = *( static_cast< int* >( mxGetData( prhs[ 0 ] ) ) );
    itk::Object::ConstPointer object = kvl::MatlabObjectArray::GetInstance()->GetObject( optimizerHandle );
    
    if ( typeid( *object ) == typeid( AtlasMeshDeformationConjugateGradientOptimizerCPU ) )
      {
      AtlasMeshDeformationConjugateGradientOptimizerCPU::ConstPointer constOptimizer 
          = static_cast< const AtlasMeshDeformationConjugateGradientOptimizerCPU* >( object.GetPointer() );
      AtlasMeshDeformationConjugateGradientOptimizerCPU::Pointer  optimizer 
          = const_cast< AtlasMeshDeformationConjugateGradientOptimizerCPU* >( constOptimizer.GetPointer() );
        
      // Set the means and variances
      optimizer->SetMeans( means );
      optimizer->SetPrecisions( precisions );
      }
    else
      {
      mexErrMsgTxt( "mesh doesn't refer to the correct ITK object type" );
      }
      
    }

protected:
  SetOptimizerPropertiesCPU() {};
  virtual ~SetOptimizerPropertiesCPU() {};


  SetOptimizerPropertiesCPU(const Self&); //purposely not implemented
  void operator=(const Self&); //purposely not implemented

private:

};

} // end namespace kvl



