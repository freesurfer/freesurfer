#include "kvlMatlabRunner.h" 
#include "kvlMatlabObjectArray.h"
#include "itkImage.h"
#include "itkImageRegionIterator.h"


namespace kvl
{

  
template< class ImageType >
class SetImageBufferHelper
{
public:

  // 
  static bool AssignImageBuffer( const mxArray*  matlabObject, itk::Object*  imageObject )
    {
    itk::Object::Pointer  itkObject = ImageConverter< ImageType >::Convert( matlabObject );
    if ( !itkObject )
      {
      return false;
      }
      
    // 
    typename ImageType::Pointer  sourceImage = static_cast< ImageType* >( itkObject.GetPointer() );
    typename ImageType::Pointer  destinationImage = static_cast< ImageType* >( imageObject );
    destinationImage->SetPixelContainer( sourceImage->GetPixelContainer() );
 
    return true;
    }
    
};
  
  
class SetImageBuffer : public MatlabRunner
{
public:
  /** Smart pointer typedef support. */
  typedef SetImageBuffer         Self;
  typedef itk::Object              Superclass;
  typedef itk::SmartPointer<Self>  Pointer;
  typedef itk::SmartPointer<const Self>  ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro( Self );

  /** Run-time type information (and related methods). */
  itkTypeMacro( SetImageBuffer, itk::Object );

  virtual void Run( int nlhs, mxArray* plhs[],
                    int nrhs, const mxArray* prhs[] )
    {
    //std::cout << "I am " << this->GetNameOfClass() 
    //          << " and I'm running! " << std::endl;
              
              
    // kvlSetImageBuffer( image, imageBuffer )
  
    // Make sure input arguments are correct
    if ( ( nrhs != 2 ) || !mxIsInt64( prhs[ 0 ] ) || 
         ( nlhs != 0 ) )
      {
      mexErrMsgTxt( "Incorrect arguments" );
      }
      
      
    // Get the object the image is pointing to
    const int imageHandle = *( static_cast< int* >( mxGetData( prhs[ 0 ] ) ) );
    itk::Object*  imageObject = kvl::MatlabObjectArray::GetInstance()->GetObject( imageHandle );

      
    // Create a new image 
    bool success = false;
    if ( !success )
      {
      success = SetImageBufferHelper< itk::Image< unsigned char, 3 > >::AssignImageBuffer( prhs[ 1 ], imageObject );
      }
    if ( !success )
      {
      success = SetImageBufferHelper< itk::Image< unsigned short, 3 > >::AssignImageBuffer( prhs[ 1 ], imageObject );
      }
    if ( !success )
      {
      success = SetImageBufferHelper< itk::Image< short, 3 > >::AssignImageBuffer( prhs[ 1 ], imageObject );
      }
    if ( !success )
      {
      success = SetImageBufferHelper< itk::Image< float, 3 > >::AssignImageBuffer( prhs[ 1 ], imageObject );
      }
    if ( !success )
      {
      std::ostringstream  errorStream;
      errorStream << "Unsupported pixel type: " << mxGetClassName( prhs[ 1 ] ); 
      mexErrMsgTxt( errorStream.str().c_str() );
      }
      

    }
  
protected:
  SetImageBuffer() {};
  virtual ~SetImageBuffer() {};


  SetImageBuffer(const Self&); //purposely not implemented
  void operator=(const Self&); //purposely not implemented

private:

};

} // end namespace kvl


