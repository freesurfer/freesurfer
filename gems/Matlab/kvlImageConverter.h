#ifndef __kvlImageConverter_h
#define __kvlImageConverter_h

#include "kvlMatlabRunner.h" 
#include "kvlMatlabObjectArray.h"
#include "itkImage.h"
#include "itkImageRegionIterator.h"
#include "itkImage.h"
#include "mex.h" 


namespace kvl
{

  
// get mxClassID of any type T
template <class T> inline mxClassID mxT() {return mxUNKNOWN_CLASS;}
template <> inline mxClassID mxT<double>() {return mxDOUBLE_CLASS;}
template <> inline mxClassID mxT<float>() {return mxSINGLE_CLASS;}
template <> inline mxClassID mxT<short>() {return mxINT16_CLASS;}
template <> inline mxClassID mxT<unsigned short>() {return mxUINT16_CLASS;}
template <> inline mxClassID mxT<unsigned char>() {return mxUINT8_CLASS;}

template< class ImageType >
class ImageConverter
{
public:

  // Convert from Matlab matrix into ITK image
  static itk::Object::Pointer  Convert( const mxArray*  matlabObject )
    {
    // Check if we can handle this matlab type
    typedef typename ImageType::PixelType  PixelType;
    mxClassID  mxclass = mxT< PixelType >();
    if ( mxGetClassID( matlabObject ) != mxclass )
      {
      return 0;  
      }
      
    // Determine the size of the image to be created  
    typedef typename ImageType::SizeType  SizeType;
    SizeType  imageSize;
    for ( int i = 0; i < 3; i++ )
      {
      imageSize[ i ] = mxGetDimensions( matlabObject )[ i ];
      //std::cout << "imageSize[ i ]: " << imageSize[ i ] << std::endl;
      }
      
      
    // Construct an ITK image
    typename ImageType::Pointer  image = ImageType::New();
    image->SetRegions( imageSize );
    image->Allocate();

    // Loop over all voxels and copy contents
    PixelType*  data = static_cast< PixelType* >( mxGetData( matlabObject ) ); 
    itk::ImageRegionIterator< ImageType >  it( image,
                                               image->GetBufferedRegion() );
    for ( ;!it.IsAtEnd(); ++it, ++data )
      {
      it.Value() = *data;
      }
    
    return image.GetPointer();  
    }  
  

  // Convert from ITK image into Matlab matrix
  static mxArray*  Convert( const itk::Object*  itkObject )
    {
    // Check if we can handle this object
    // if ( typeid( *itkObject ) != typeid( ImageType ) )
     if ( strcmp(typeid( *itkObject ).name(), typeid( ImageType ).name()) )  // Eugenio: MAC compatibility
      {
      return 0;
      }

    typedef typename ImageType::PixelType  PixelType;

    // Get the image
    typename ImageType::ConstPointer  image 
          = static_cast< const ImageType* >( itkObject );
          
    // Get it's size
    mwSize  dims[ 3 ];
    for ( int i = 0; i < 3; i++ )
      {
      dims[ i ] = image->GetBufferedRegion().GetSize()[ i ];
      }
      
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
  
    return matlabObject;  
    }


};

  
} // end namespace kvl

#endif
