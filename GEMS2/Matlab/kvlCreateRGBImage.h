#include "kvlMatlabRunner.h" 
#include "kvlMatlabObjectArray.h"
#include "itkImage.h"
#include "itkImageRegionIterator.h"
#include "itkRGBAPixel.h"


namespace kvl
{

  
class CreateRGBImage : public MatlabRunner
{
public:
  /** Smart pointer typedef support. */
  typedef CreateRGBImage         Self;
  typedef itk::Object              Superclass;
  typedef itk::SmartPointer<Self>  Pointer;
  typedef itk::SmartPointer<const Self>  ConstPointer;
  typedef itk::RGBAPixel< unsigned char > RGBAPixelType;

  /** Method for creation through the object factory. */
  itkNewMacro( Self );

  /** Run-time type information (and related methods). */
  itkTypeMacro( CreateRGBImage, itk::Object );

  virtual void Run( int nlhs, mxArray* plhs[],
                    int nrhs, const mxArray* prhs[] )
    {
    //std::cout << "I am " << this->GetNameOfClass() 
    //          << " and I'm running! " << std::endl;
              
              
    // image = kvlCreateImage( r,g,b,a )
  
    // Make sure input arguments are correct
    if ( ( nrhs != 4 ) || ( nlhs != 1 ) )
      {
      mexErrMsgTxt( "Incorrect arguments" );
      }

    if(!mxIsUint8(prhs[0]) || !mxIsUint8(prhs[1]) || !mxIsUint8(prhs[2]) || !mxIsUint8(prhs[3])){
       mexErrMsgTxt( "Unsupported data type" );
    }
      
   
     // Determine the size of the image to be created
    typedef itk::Image< RGBAPixelType, 3 >   RGBAImageType;   
    typedef RGBAImageType::SizeType  SizeType;
    SizeType  imageSize;
    for ( int i = 0; i < 3; i++ )
      {
      imageSize[ i ] = mxGetDimensions( prhs[ 0 ] )[ i ];
      std::cout << "imageSize[ i ]: " << imageSize[ i ] << std::endl;
      }
      
      
    // Construct an ITK image
    RGBAImageType::Pointer  RGBAImage = RGBAImageType::New();
    RGBAImage->SetRegions( imageSize );
    RGBAImage->Allocate();
 
    // Loop over all voxels and copy contents
    

    unsigned char*  data1 = static_cast< unsigned char* >( mxGetData( prhs[0] ) );
    unsigned char*  data2 = static_cast< unsigned char* >( mxGetData( prhs[1] ) );
    unsigned char*  data3 = static_cast< unsigned char* >( mxGetData( prhs[2] ) );
    unsigned char*  data4 = static_cast< unsigned char* >( mxGetData( prhs[3] ) );  
    itk::ImageRegionIterator< RGBAImageType >  it( RGBAImage,
                                               RGBAImage->GetBufferedRegion() );
    RGBAPixelType pixel;
    int counter = 1;
    for ( ;!it.IsAtEnd(); ++it, ++data1, ++data2, ++data3, ++data4 )
      {

      if(counter % 10000 == 0){
         std::cout << "pixelValues: " << static_cast<unsigned>(*data1) <<","<< static_cast<unsigned>(*data2) << "," << static_cast<unsigned>(*data3) << "," << static_cast<unsigned>(*data4)<< std::endl;
      }
      pixel.SetRed(*data1);
      pixel.SetGreen(*data2);
      pixel.SetBlue(*data3);
      pixel.SetAlpha(*data4);
      it.Set(pixel);
      counter += 1;
      }
    
  

    
    // Store the created image in persistent memory
    const int  imageHandle = kvl::MatlabObjectArray::GetInstance()->AddObject( RGBAImage.GetPointer() );
     
    // Return the handle to Matlab
    mwSize  dims[ 1 ];
    dims[ 0 ] = 1;
    plhs[ 0 ] = mxCreateNumericArray( 1, dims, mxINT64_CLASS, mxREAL );
    *( static_cast< int* >( mxGetData( plhs[ 0 ] ) ) ) = imageHandle;


    }  
  
protected:
  CreateRGBImage() {};
  virtual ~CreateRGBImage() {};


  CreateRGBImage(const Self&); //purposely not implemented
  void operator=(const Self&); //purposely not implemented

private:

};

} // end namespace kvl


