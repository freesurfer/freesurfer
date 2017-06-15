#include "kvlMatlabRunner.h" 
#include "kvlMatlabObjectArray.h"
#include "kvlCroppedImageReader.h"
#include "itkCastImageFilter.h"


namespace kvl
{


class GetCroppedRegion : public MatlabRunner
{
public:
  /** Smart pointer typedef support. */
  typedef GetCroppedRegion         Self;
  typedef itk::Object              Superclass;
  typedef itk::SmartPointer<Self>  Pointer;
  typedef itk::SmartPointer<const Self>  ConstPointer;
  typedef itk::Image< float, 3 >  ImageType;

  /** Method for creation through the object factory. */
  itkNewMacro( Self );

  /** Run-time type information (and related methods). */
  itkTypeMacro( GetCroppedRegion, itk::Object );

  virtual void Run( int nlhs, mxArray* plhs[],
                    int nrhs, const mxArray* prhs[] )
    {
    //std::cout << "I am " << this->GetNameOfClass() 
    //          << " and I'm running! " << std::endl;
              
    // [origIndex croppedIndex croppedSize origSize] = kvlGetCroppedRange( imageFileName, boundingFileName )
  
    // Make sure input arguments are correct
    if ( ( nrhs != 2 ) || 
         !mxIsChar( prhs[ 0 ] ) || 
         !mxIsChar( prhs[ 1 ] ) || 
         ( nlhs != 4 ) )
      {
      mexErrMsgTxt( "Incorrect arguments" );
      }

    // Retrieve input arguments
    const std::string imageFileName = mxArrayToString( prhs[0] );
    const std::string boundingFileName = mxArrayToString( prhs[1] );

    // Read the image
    kvl::CroppedImageReader::Pointer  reader = kvl::CroppedImageReader::New();
    //reader->SetExtraFraction( 0.1 );
    reader->SetExtraFraction( 0.0 );
    reader->Read( imageFileName.c_str(), boundingFileName.c_str() );
 
    //Get the cropped region in the original image
    ImageType::RegionType  originalImageOriginalRegion;    
    originalImageOriginalRegion = reader->GetOriginalImageOriginalRegion();
    //ImageType::IndexType origRegionIndex;
    //origRegionIndex = originalImageOriginalRegion.GetIndex();

    ImageType::SizeType origArea;
    origArea = originalImageOriginalRegion.GetSize();

    ImageType::RegionType  originalImageCroppedRegion;    
    originalImageCroppedRegion = reader->GetOriginalImageRegion();
    ImageType::IndexType croppedRegionIndex;
    croppedRegionIndex = originalImageCroppedRegion.GetIndex();

    ImageType::SizeType croppedArea;
    croppedArea = originalImageCroppedRegion.GetSize();

    ImageType::RegionType  CroppedRegion;
    CroppedRegion = reader->GetCroppedImageRegion();
    ImageType::IndexType croppedIndex;
    croppedIndex = CroppedRegion.GetIndex();
    

    
    // Return the data to Matlab

    // First the cropping indices
    mwSize  dims[ 2 ];
    dims[ 0 ] = 1;
    dims[ 1 ] = 3;
    plhs[ 0 ] = mxCreateNumericArray( 2, dims, mxINT32_CLASS, mxREAL );
    int* data = static_cast< int* >( mxGetData( plhs[ 0 ]) );
    *(data) = croppedRegionIndex[ 0 ];
    *(data + 1) = croppedRegionIndex[ 1 ];
    *(data + 2) = croppedRegionIndex[ 2 ];
    
    plhs[ 1 ] = mxCreateNumericArray( 2, dims, mxINT32_CLASS, mxREAL );
    data = static_cast< int* >( mxGetData( plhs[ 1 ]) );
    *(data) = croppedIndex[ 0 ];
    *(data + 1) = croppedIndex[ 1 ];
    *(data + 2) = croppedIndex[ 2 ];

    // Then the cropped size
    plhs[ 2 ] = mxCreateNumericArray( 2, dims, mxINT32_CLASS, mxREAL );
    data = static_cast< int* >( mxGetData( plhs[ 2 ]) );
    *(data) = croppedArea[ 0 ];
    *(data + 1) = croppedArea[ 1 ];
    *(data + 2) = croppedArea[ 2 ];

    // Finally the original size
    plhs[ 3 ] = mxCreateNumericArray( 2, dims, mxINT32_CLASS, mxREAL );
    data = static_cast< int* >( mxGetData( plhs[ 3 ]) );
    *(data) = origArea[ 0 ];
    *(data + 1) = origArea[ 1 ];
    *(data + 2) = origArea[ 2 ];

    }
  
protected:
  GetCroppedRegion() {};
  virtual ~GetCroppedRegion() {};


  GetCroppedRegion(const Self&); //purposely not implemented
  void operator=(const Self&); //purposely not implemented

private:

};

} // end namespace kvl

