#include "kvlMatlabRunner.h" 
#include "kvlMatlabObjectArray.h"
#include "kvlCroppedImageReader.h"
#include "itkCastImageFilter.h"


namespace kvl
{


class ReadImage : public MatlabRunner
{
public:
  /** Smart pointer typedef support. */
  typedef ReadImage         Self;
  typedef itk::Object              Superclass;
  typedef itk::SmartPointer<Self>  Pointer;
  typedef itk::SmartPointer<const Self>  ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro( Self );

  /** Run-time type information (and related methods). */
  itkTypeMacro( ReadImage, itk::Object );

  virtual void Run( int nlhs, mxArray* plhs[],
                    int nrhs, const mxArray* prhs[] )
    {
    //std::cout << "I am " << this->GetNameOfClass() 
    //          << " and I'm running! " << std::endl;
              
    // [ image, transform ] = kvlReadImage( imageFileName )
  
    // Make sure input arguments are correct
    if ( ( nrhs != 1 ) || 
         !mxIsChar( prhs[ 0 ] ) || 
         ( nlhs != 2 ) )
      {
      mexErrMsgTxt( "Incorrect arguments" );
      }

    // Retrieve input arguments
    const std::string imageFileName = mxArrayToString( prhs[0] );

    // Read the image
    kvl::CroppedImageReader::Pointer  reader = kvl::CroppedImageReader::New();
    reader->Read( imageFileName.c_str() );

    // Convert the image to float
    typedef itk::Image< float, 3 >  ImageType;
    typedef itk::CastImageFilter< kvl::CroppedImageReader::ImageType, ImageType >  CasterType;
    CasterType::Pointer  caster = CasterType::New();
    caster->SetInput( reader->GetImage() );
    caster->Update();

    
    // Store the image and transform in persistent memory
    const int imageHandle = kvl::MatlabObjectArray::GetInstance()->AddObject( caster->GetOutput() );
    typedef kvl::CroppedImageReader::TransformType  TransformType;
    TransformType::Pointer  transform = TransformType::New();
    reader->GetWorldToImageTransform()->GetInverse( transform );
    const int transformHandle = kvl::MatlabObjectArray::GetInstance()->AddObject( transform );
    
    // Return the handles to Matlab
    mwSize  dims[ 1 ];
    dims[ 0 ] = 1;
    plhs[ 0 ] = mxCreateNumericArray( 1, dims, mxINT64_CLASS, mxREAL );
    *( static_cast< int* >( mxGetData( plhs[ 0 ] ) ) ) = imageHandle;
    plhs[ 1 ] = mxCreateNumericArray( 1, dims, mxINT64_CLASS, mxREAL );
    *( static_cast< int* >( mxGetData( plhs[ 1 ] ) ) ) = transformHandle;
    }
  
protected:
  ReadImage() {};
  virtual ~ReadImage() {};


  ReadImage(const Self&); //purposely not implemented
  void operator=(const Self&); //purposely not implemented

private:

};

} // end namespace kvl

