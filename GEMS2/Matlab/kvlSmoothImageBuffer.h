#include "kvlMatlabRunner.h" 
#include "kvlMatlabObjectArray.h"
#include "itkImage.h"
#include "itkImageRegionIterator.h"
#include "itkImageRegionConstIterator.h"
#include "itkDiscreteGaussianImageFilter.h"
#include "itkCastImageFilter.h"
#include "kvlImageConverter.h"



namespace kvl
{
  
  
template< class ImageType >
class SmoothImageBufferHelper
{
public:

  // 
  static mxArray* Smooth( const mxArray*  matlabObject, const double* sigmas )
    {
    // Check if we can handle this type
    ImageConverter< ImageType >  converter;
    itk::Object::Pointer  object = converter.Convert( matlabObject );
    if ( !object )
      {
      return 0;  
      }
    
    // Get the ITK image
    typename ImageType::Pointer  image = static_cast< ImageType* >( object.GetPointer() );
      
    // Do the actual work in ITK
    typedef itk::Image< float, 3 >  InternalImageType;
    typedef itk::CastImageFilter< ImageType, InternalImageType >   CasterType;
    typedef itk::DiscreteGaussianImageFilter< InternalImageType, InternalImageType >  SmootherType;
    typedef itk::CastImageFilter< InternalImageType, ImageType >  BackCasterType;

    typename CasterType::Pointer caster = CasterType::New();
    caster->SetInput( image );
    SmootherType::Pointer smoother = SmootherType::New();
    smoother->SetInput( caster->GetOutput() );
    smoother->SetMaximumError( 0.1 );
    smoother->SetUseImageSpacingOff();
    double  variances[ 3 ];
    for ( int i = 0; i < 3; i++ )
      {
      variances[ i ] = sigmas[ i ] * sigmas[ i ];  
      }  
    smoother->SetVariance( variances );
    typename BackCasterType::Pointer  backCaster = BackCasterType::New();
    backCaster->SetInput( smoother->GetOutput() );
    backCaster->Update();
    typename ImageType::ConstPointer  smoothedImage = backCaster->GetOutput();
 
    // Convert back to Matlab
    return converter.Convert( smoothedImage );
    }
  
};

  

class SmoothImageBuffer : public MatlabRunner
{
public:
  /** Smart pointer typedef support. */
  typedef SmoothImageBuffer         Self;
  typedef itk::Object              Superclass;
  typedef itk::SmartPointer<Self>  Pointer;
  typedef itk::SmartPointer<const Self>  ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro( Self );

  /** Run-time type information (and related methods). */
  itkTypeMacro( SmoothImageBuffer, itk::Object );

  virtual void Run( int nlhs, mxArray* plhs[],
                    int nrhs, const mxArray* prhs[] )
    {
    //std::cout << "I am " << this->GetNameOfClass() 
    //          << " and I'm running! " << std::endl;
              
              
    // smoothedImageBuffer = kvlSmoothImageBuffer( imageBuffer, sigma(s) )
 
    // Make sure input arguments are correct
    if ( ( nrhs != 2 ) || ( nlhs != 1 ) || !mxIsDouble( prhs[ 1 ] ) )
      {
      mexErrMsgTxt( "Incorrect arguments" );
      }
      
    // Retrieve dimensions of the input image
    if ( mxGetNumberOfDimensions( prhs[ 0 ] ) != 3 )
      {
      mexErrMsgTxt( "Input must be 3-dimensional real matrix" );
      }
     
    // Retrieve the smoothing sigmas
    double  sigmas[ 3 ];
    double*  tmp = mxGetPr( prhs[ 1 ] );
    const int  numberOfRows = *( mxGetDimensions( prhs[ 1 ] ) );
    const int  numberOfColumns = *( mxGetDimensions( prhs[ 1 ] ) + 1 );
    //std::cout << "numberOfRows: " << numberOfRows << std::endl;
    //std::cout << "numberOfColumns: " << numberOfColumns << std::endl;
    if ( ( numberOfRows * numberOfColumns ) == 1 )
      {
      for ( int i = 0; i < 3; i++ )
        {
        sigmas[ i ] = *tmp;
        }
        
      }
    else if ( ( numberOfRows * numberOfColumns ) == 3 )
      {
      for ( int i = 0; i < 3; i++, tmp++ )
        {
        sigmas[ i ] = *tmp;
        }
      }
    else
      {
      mexErrMsgTxt( "Sigma(s) must be 1- or 3-dimensional" );
      }
    
    
    
    mxArray*  matlabObject = 0;
    if ( !matlabObject )
      {
      matlabObject = SmoothImageBufferHelper< itk::Image< unsigned char, 3 > >::Smooth( prhs[ 0 ], sigmas );
      }
    if ( !matlabObject )
      {
      matlabObject = SmoothImageBufferHelper< itk::Image< unsigned short, 3 > >::Smooth( prhs[ 0 ], sigmas );
      }
    if ( !matlabObject )
      {
      matlabObject = SmoothImageBufferHelper< itk::Image< short, 3 > >::Smooth( prhs[ 0 ], sigmas );
      }
    if ( !matlabObject )
      {
      matlabObject = SmoothImageBufferHelper< itk::Image< float, 3 > >::Smooth( prhs[ 0 ], sigmas );
      }
    if ( !matlabObject )
      {
      mexErrMsgTxt( "Unsupported pixel type" );
      }

    plhs[ 0 ] = matlabObject; 
    }
  
protected:
  SmoothImageBuffer() {};
  virtual ~SmoothImageBuffer() {};


  SmoothImageBuffer(const Self&); //purposely not implemented
  void operator=(const Self&); //purposely not implemented

private:

};

} // end namespace kvl







