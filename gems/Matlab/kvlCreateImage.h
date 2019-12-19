#include "kvlMatlabRunner.h" 
#include "kvlMatlabObjectArray.h"
#include "kvlImageConverter.h"
#include "itkImage.h"
#include "itkImageRegionIterator.h"
#include "itkRGBPixel.h"


namespace kvl
{

  
class CreateImage : public MatlabRunner
{
public:
  /** Smart pointer typedef support. */
  typedef CreateImage         Self;
  typedef itk::Object              Superclass;
  typedef itk::SmartPointer<Self>  Pointer;
  typedef itk::SmartPointer<const Self>  ConstPointer;
  typedef itk::RGBPixel< float > RGBPixelType;

  /** Method for creation through the object factory. */
  itkNewMacro( Self );

  /** Run-time type information (and related methods). */
  itkTypeMacro( CreateImage, itk::Object );

  virtual void Run( int nlhs, mxArray* plhs[],
                    int nrhs, const mxArray* prhs[] )
    {
    //std::cout << "I am " << this->GetNameOfClass() 
    //          << " and I'm running! " << std::endl;
              
              
    // image = kvlCreateImage( imageBuffer )
  
    // Make sure input arguments are correct
    if ( ( nrhs != 1 ) || ( nlhs != 1 ) )
      {
      mexErrMsgTxt( "Incorrect arguments" );
      }
      
    // Retrieve dimensions of the input image
    //if ( mxGetNumberOfDimensions( prhs[ 0 ] ) != 3 )
    //  {
    //  mexErrMsgTxt( "Input must be 3-dimensional real matrix" );
    //  }
     
     
      
    //
    itk::Object::Pointer  itkObject = 0;

    if( mxGetNumberOfDimensions( prhs[ 0 ] ) == 3 ){
    	if ( !itkObject )
      	{
      	ImageConverter< itk::Image< unsigned char, 3 > >  converter;
      	itkObject = converter.Convert( prhs[ 0 ] );
      	}
    	if ( !itkObject )
      	{
      	ImageConverter< itk::Image< unsigned short, 3 > >  converter;
      	itkObject = converter.Convert( prhs[ 0 ] );
      	}
    	if ( !itkObject )
      	{
      	ImageConverter< itk::Image< short, 3 > >  converter;
      	itkObject = converter.Convert( prhs[ 0 ] );
      	}
    	if ( !itkObject )
     	{
      	ImageConverter< itk::Image< float, 3 > >  converter;
      	itkObject = converter.Convert( prhs[ 0 ] );
      	}
    	if ( !itkObject )
     	{
      	std::ostringstream  errorStream;
      	errorStream << "Unsupported pixel type: " << mxGetClassName( prhs[ 0 ] ); 
      	mexErrMsgTxt( errorStream.str().c_str() );
      	}
    }
    else if( mxGetNumberOfDimensions( prhs[ 0 ] ) == 4){
	if ( !itkObject )
      	{
      	ImageConverter< itk::Image< RGBPixelType, 3 > >  converter;
      	itkObject = converter.Convert( prhs[ 0 ] );
        std::cout<<"Creating a color image!"<<std::endl;
      	}
    }




    
    // Store the created image in persistent memory
    const int  imageHandle = kvl::MatlabObjectArray::GetInstance()->AddObject( itkObject );
     
    // Return the handle to Matlab
    mwSize  dims[ 1 ];
    dims[ 0 ] = 1;
    plhs[ 0 ] = mxCreateNumericArray( 1, dims, mxINT64_CLASS, mxREAL );
    *( static_cast< int* >( mxGetData( plhs[ 0 ] ) ) ) = imageHandle;


    }  
  
protected:
  CreateImage() {};
  virtual ~CreateImage() {};


  CreateImage(const Self&); //purposely not implemented
  void operator=(const Self&); //purposely not implemented

private:

};

} // end namespace kvl


