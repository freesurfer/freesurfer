#include "kvlMatlabRunner.h" 
#include "kvlMatlabObjectArray.h"
#include "kvlImageConverter.h"
#include "itkImage.h"
#include "itkImageRegionConstIterator.h"


namespace kvl
{

  
class GetImageBuffer : public MatlabRunner
{
public:
  /** Smart pointer typedef support. */
  typedef GetImageBuffer         Self;
  typedef itk::Object              Superclass;
  typedef itk::SmartPointer<Self>  Pointer;
  typedef itk::SmartPointer<const Self>  ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro( Self );

  /** Run-time type information (and related methods). */
  itkTypeMacro( GetImageBuffer, itk::Object );

  virtual void Run( int nlhs, mxArray* plhs[],
                    int nrhs, const mxArray* prhs[] )
    {
    //std::cout << "I am " << this->GetNameOfClass() 
    //          << " and I'm running! " << std::endl;
              
              
      // imageBuffer = kvlGetImageBuffer( image )
  
    // Make sure input arguments are correct
    if ( ( nrhs != 1 ) || !mxIsInt64( prhs[ 0 ] ) || 
         ( nlhs != 1 ) )
      {
      mexErrMsgTxt( "Incorrect arguments" );
      }
      
    // Retrieve input arguments
    const int imageHandle = *( static_cast< int* >( mxGetData( prhs[ 0 ] ) ) );

    // Get the object
    itk::Object::ConstPointer object = kvl::MatlabObjectArray::GetInstance()->GetObject( imageHandle );

    // Copy the buffer contents into a Matlab matrix
    mwSize  dims[ 1 ];
    dims[ 0 ] = 1;

    mxArray*  matlabObject = 0;
    if ( !matlabObject )
      {
      ImageConverter< itk::Image< unsigned char, 3 > >  converter;
      matlabObject = converter.Convert( object.GetPointer() );
      }
    if ( !matlabObject )
      {
      ImageConverter< itk::Image< unsigned short, 3 > >  converter;
      matlabObject = converter.Convert( object.GetPointer() );
      }
    if ( !matlabObject )
      {
      ImageConverter< itk::Image< short, 3 > >  converter;
      matlabObject = converter.Convert( object.GetPointer() );
      }
    if ( !matlabObject )
      {
      ImageConverter< itk::Image< float, 3 > >  converter;
      matlabObject = converter.Convert( object.GetPointer() );
      }
    if ( !matlabObject )
      {
      mexErrMsgTxt( "Unsupported pixel type" );
      }

    plhs[ 0 ] = matlabObject; 

    }
  
protected:
  GetImageBuffer() {};
  virtual ~GetImageBuffer() {};


  GetImageBuffer(const Self&); //purposely not implemented
  void operator=(const Self&); //purposely not implemented

private:
  


};

} // end namespace kvl


