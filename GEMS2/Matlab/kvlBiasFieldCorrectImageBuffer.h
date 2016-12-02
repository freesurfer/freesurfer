#include "kvlMatlabRunner.h" 
#include "kvlMatlabObjectArray.h"
#include "itkImage.h"
#include "itkImageRegionIterator.h"
#include "itkImageRegionConstIterator.h"
#include "kvlEMSegmenter.h"


namespace kvl
{

class BiasFieldCorrectImageBuffer : public MatlabRunner
{
public:
  /** Smart pointer typedef support. */
  typedef BiasFieldCorrectImageBuffer         Self;
  typedef itk::Object              Superclass;
  typedef itk::SmartPointer<Self>  Pointer;
  typedef itk::SmartPointer<const Self>  ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro( Self );

  /** Run-time type information (and related methods). */
  itkTypeMacro( BiasFieldCorrectImageBuffer, itk::Object );

  virtual void Run( int nlhs, mxArray* plhs[],
                    int nrhs, const mxArray* prhs[] )
    {
    std::cout << "I am " << this->GetNameOfClass() 
              << " and I'm running! " << std::endl;
              
              
    // [ biasFieldCorrectedImageBuffer, biasFieldImageBuffer ] = kvlBiasFieldCorrectImageBuffer( imageBuffer, biasFieldOrder )
 
    // Some typedefs
    typedef itk::Image< unsigned short, 3 >   ImageType;

              
    // Make sure input arguments are correct
    if ( ( nrhs != 2 ) || ( mxGetClassID( prhs[ 0 ] ) != mxUINT16_CLASS ) || !mxIsDouble( prhs[ 1 ] ) )
      {
      mexErrMsgTxt( "Incorrect arguments" );
      }
      
    // Retrieve dimensions of the input image
    if ( mxGetNumberOfDimensions( prhs[ 0 ] ) != 3 )
      {
      mexErrMsgTxt( "Image buffer must be 3-dimensional real matrix" );
      }
     
    typedef FImageType::SizeType  SizeType;
    SizeType  imageSize;
    imageSize[ 0 ] = mxGetDimensions( prhs[ 0 ] )[ 0 ];
    imageSize[ 1 ] = mxGetDimensions( prhs[ 0 ] )[ 1 ];
    imageSize[ 2 ] = mxGetDimensions( prhs[ 0 ] )[ 2 ];
    

    // Retrieve the smoothing sigma
    const int biasFieldOrder = static_cast< int >( *( mxGetPr( prhs[ 1 ] ) ) );
    std::cout << "imageSize: " << imageSize << std::endl;
    std::cout << "biasFieldOrder: " << biasFieldOrder << std::endl;
    
        
    // Construct an ITK image
    ImageType::Pointer  image = ImageType::New();
    image->SetRegions( imageSize );
    image->Allocate();
    ImageType::PixelType*  data 
      = static_cast< ImageType::PixelType* >( mxGetData( prhs[ 0 ] ) ); 
    itk::ImageRegionIterator< ImageType >  it( image,
                                               image->GetBufferedRegion() );
    for ( ;!it.IsAtEnd(); ++it, ++data )
      {
      it.Value() = *data;
      }
        
        
    // Do the actual work in ITK



        // Create the correct Matlab matrix type 
        mwSize  dims[ 3 ];
        for ( int i = 0; i < 3; i++ )
          {
          dims[ i ] = imageSize[ i ];
          }        
        plhs[ 0 ] = mxCreateNumericArray( 3, dims, mxSINGLE_CLASS, mxREAL );
        data = static_cast< ImageType::PixelType* >( mxGetData( plhs[ 0 ] ) ); 

        // Loop over all voxels and copy contents
        itk::ImageRegionConstIterator< ImageType >  it2( smoothedImage,
                                                         smoothedImage->GetBufferedRegion() );
        for ( ;!it2.IsAtEnd(); ++it2, ++data )
          {
          *data = it2.Value();
          }
    
        }
        break;
      default:
        std::ostringstream  errorStream;
        errorStream << "Unsupported pixel type: " << mxGetClassName( prhs[ 0 ] ); 
        mexErrMsgTxt( errorStream.str().c_str() );
      } // End check pixel type


    }
  
protected:
  BiasFieldCorrectImageBuffer() {};
  virtual ~BiasFieldCorrectImageBuffer() {};


  BiasFieldCorrectImageBuffer(const Self&); //purposely not implemented
  void operator=(const Self&); //purposely not implemented

private:

};

} // end namespace kvl







