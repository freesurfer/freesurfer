#include "kvlMatlabRunner.h" 
#include "kvlMatlabObjectArray.h"
#include "kvlImageConverter.h"
#include "itkImage.h"
#include "itkImageRegionConstIterator.h"



namespace kvl
{

  
class GetImageFromOptimizer : public MatlabRunner
{
public:
  /** Smart pointer typedef support. */
  typedef GetImageFromOptimizer         Self;
  typedef itk::Object              Superclass;
  typedef itk::SmartPointer<Self>  Pointer;
  typedef itk::SmartPointer<const Self>  ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro( Self );

  /** Run-time type information (and related methods). */
  itkTypeMacro( GetImageFromOptimizer, itk::Object );

  virtual void Run( int nlhs, mxArray* plhs[],
                    int nrhs, const mxArray* prhs[] )
    {
    //std::cout << "I am " << this->GetNameOfClass() 
    //          << " and I'm running! " << std::endl;
              
              
    // imageBuffer = kvlGetImageFromOptimizer( optimizer, ind )
  
    // Make sure input arguments are correct
    if ( ( nrhs < 2 ) || 
         !mxIsInt64( prhs[ 0 ] ))// || 
         //( nlhs != 1 ) )
      {
      mexErrMsgTxt( "Incorrect arguments" );
      }
      
        
    typedef itk::Image< float, 3 >  ImageType;
    // Retrieve the optimizer, and set the means and variances
    const int optimizerHandle = *( static_cast< int* >( mxGetData( prhs[ 0 ] ) ) );
    itk::Object::ConstPointer object = kvl::MatlabObjectArray::GetInstance()->GetObject( optimizerHandle );
    const int imageNumber = mxGetScalar(prhs[1]);
    itk::Image< float, 3 >::ConstPointer image;
    
    //std::cout<<"TypeId of object: "<<typeid(*object)<<std::endl;
    //std::cout<<"TypeId of LM: "<<typeid(AtlasMeshDeformationLevenbergMarquardtOptimizer)<<std::endl;
    //std::cout<<"TypeId of CG: "<<typeid(AtlasMeshDeformationConjugateGradientOptimizer)<<std::endl;
    //std::cout<<"TypeId of CGMulti: "<<typeid(AtlasMeshDeformationConjugateGradientOptimizerMultiAtlas)<<std::endl;
    if ( typeid( *object ) == typeid( AtlasMeshDeformationLevenbergMarquardtOptimizer ) )
      {
      AtlasMeshDeformationLevenbergMarquardtOptimizer::ConstPointer constOptimizer 
          = static_cast< const AtlasMeshDeformationLevenbergMarquardtOptimizer* >( object.GetPointer() );
      AtlasMeshDeformationLevenbergMarquardtOptimizer::Pointer  optimizer 
          = const_cast< AtlasMeshDeformationLevenbergMarquardtOptimizer* >( constOptimizer.GetPointer() );
        
      // Get the image from optimizer
      //itk::Object::ConstPointer imageObject = optimizer->GetImage(imageNumber);
      image = optimizer->GetImage(imageNumber);
      }
    else if ( typeid( *object ) == typeid( AtlasMeshDeformationConjugateGradientOptimizer ) )
      {
      AtlasMeshDeformationConjugateGradientOptimizer::ConstPointer constOptimizer 
          = static_cast< const AtlasMeshDeformationConjugateGradientOptimizer* >( object.GetPointer() );
      AtlasMeshDeformationConjugateGradientOptimizer::Pointer  optimizer 
          = const_cast< AtlasMeshDeformationConjugateGradientOptimizer* >( constOptimizer.GetPointer() );
        
      // Get the image from optimizer
      //itk::Object::ConstPointer imageObject = optimizer->GetImage(imageNumber);
      image = optimizer->GetImage(imageNumber);
      }
    else
      {
      mexErrMsgTxt( "mesh doesn't refer to the correct ITK object type" );
      }

    // Get it's size
    mwSize  dims[ 3 ];
    for ( int i = 0; i < 3; i++ )
      {
      //dims[ i ] = image->GetBufferedRegion().GetSize()[ i ];
      dims[ i ] = image->GetLargestPossibleRegion().GetSize()[ i ];
      }
      

    typedef ImageType::PixelType  PixelType;
    // Create the correct Matlab matrix type
    mxClassID mxclass = mxT< PixelType >();
    mxArray* matlabObject = mxCreateNumericArray( 3, dims, mxclass, mxREAL );
    PixelType*  data 
        = static_cast< PixelType* >( mxGetData( matlabObject ) );  

    // Loop over all voxels and copy contents
    itk::ImageRegionConstIterator< ImageType >  it( image,
                                                    image->GetBufferedRegion() );
    for ( ;!it.IsAtEnd(); ++it, ++data )
      {
      *data = it.Value();
      }

    mexPrintf("dimensions = %d\n",dims[0]);
    mexPrintf("dimensions = %d\n",dims[1]);
    mexPrintf("dimensions = %d\n",dims[2]);
    plhs[ 0 ] = matlabObject;
      
    }

protected:
  GetImageFromOptimizer() {};
  virtual ~GetImageFromOptimizer() {};


  GetImageFromOptimizer(const Self&); //purposely not implemented
  void operator=(const Self&); //purposely not implemented

private:

};

} // end namespace kvl



