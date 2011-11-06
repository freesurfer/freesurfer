#ifndef _SBL_IMAGE_TRANSFORM_H_
#define _SBL_IMAGE_TRANSFORM_H_
#include <sbl/core/Pointer.h>
#include <sbl/core/File.h>
#include <sbl/math/Vector.h>
#include <sbl/math/Geometry.h>
#include <sbl/image/Image.h>
namespace sbl {


/*! \file ImageTransform.h
    \brief The ImageTransform module provides functions for spatial transformations of images.
    Most of these functions are    simple wrappers for OpenCV functions.  
    For other image filters, see the ImageUtil module.
*/


//-------------------------------------------
// IMAGE TRANSFORMATION FUNCTIONS
//-------------------------------------------


/// extract sub-image
template <typename ImageType> aptr<ImageType> crop( const ImageType &input, int xMin, int xMax, int yMin, int yMax );


/// resize image
template <typename ImageType> aptr<ImageType> resize( const ImageType &input, int newWidth, int newHeight, bool filter );


/// translate and scale an image
template <typename ImageType> aptr<ImageType> shiftScale( const ImageType &input, float xOffset, float yOffset, float xScale, float yScale, int outputWidth, int outputHeight );


/// apply linear transformation to image
template <typename ImageType> aptr<ImageType> warpAffine( const ImageType &input, float xOffset, float yOffset, float x1, float y1, float x2, float y2, int outputWidth, int outputHeight, int fillColor );


/// flip image vertically (about horizontal axis)
template <typename ImageType> aptr<ImageType> flipVert( const ImageType &input );


/// flip image horizontally (about vertical axis)
template <typename ImageType> aptr<ImageType> flipHoriz( const ImageType &input );
//aptr<ImageGrayU> flipHoriz( const ImageGrayU &input );


/// rotate image 180 degrees
//aptr<ImageGrayU> rotate180( const ImageGrayU &input );
template <typename ImageType> aptr<ImageType> rotate180( const ImageType &input );

/// rotate 90 degrees (counter-clockwise)
//aptr<ImageGrayU> rotate90( ImageGrayU &input );
template <typename ImageType> aptr<ImageType> rotate90( const ImageType &input );


/// rotate 270 degrees (counter-clockwise)
//aptr<ImageGrayU> rotate270( ImageGrayU &input );
template <typename ImageType> aptr<ImageType> rotate270( const ImageType &input );


/// rotate arbitrary amount (counter-clockwise)
//aptr<ImageGrayU> rotate( ImageGrayU &input, float angleDegrees, int fillColor );
template <typename ImageType> aptr<ImageType> rotate( const ImageType &input, float angleDegrees, int fillColor  );


//-------------------------------------------
// IMAGE TRANSFORM CLASS
//-------------------------------------------


/// The ImageTransform class represents a 2D image transformation (currently only linear transformations).
class ImageTransform {
public:

    /// create transformation from parameters (assumes correct number of parameters)
    explicit ImageTransform( const VectorF &params ) : m_params( params ) {}

    /// create a translation transformation
    ImageTransform( float xOffset, float yOffset );

    /// create an affine transformation
    ImageTransform( float xOffset, float yOffset, float xScale, float yScale );

    /// load transformation parameters from file
    explicit ImageTransform( File &file );

    /// save transformation parameters to file
    void save( File &file ) const;

    /// transform a point
    inline float xTransform( float x, float y ) const { return (m_params.length() == 2) ? m_params[ 0 ] + x : m_params[ 0 ] + x * m_params[ 2 ] + y * m_params[ 4 ]; }
    inline float yTransform( float x, float y ) const { return (m_params.length() == 2) ? m_params[ 1 ] + y : m_params[ 1 ] + x * m_params[ 3 ] + y * m_params[ 5 ]; }

    /// get a single transformation parameter
    inline float param( int index ) const { return m_params[ index ]; }

    /// get number of parameters
    inline int paramCount() const { return m_params.length(); }

    /// get the translation components of the transformation
    inline float xOffset() const { return m_params[ 0 ]; }
    inline float yOffset() const { return m_params[ 1 ]; }

    /// set the translation components of the transformation
    inline void setOffset( float xOffset, float yOffset ) { m_params[ 0 ] = xOffset; m_params[ 1 ] = yOffset; }

    /// compute the inverse of the transformation
    aptr<ImageTransform> inverse() const;

    /// map an image forward according to the transformation
    aptr<ImageGrayU> mapForward( const ImageGrayU &img, int outputWidth, int outputHeight, int fillColor ) const;
    aptr<ImageColorU> mapForward( const ImageColorU &img, int outputWidth, int outputHeight, int fillColor ) const;

    /// map an image backward using the inverse transformation
    aptr<ImageGrayU> mapBackward( const ImageGrayU &img, int outputWidth, int outputHeight, int fillColor ) const;

    /// map a point forward according to the transformation
    Point2 mapForward( const Point2 &pt ) const;

    /// map a point backward according to the inverse transformation
    Point2 mapBackward( const Point2 &pt ) const;

    /// display transformation parameters
    void display( int indent );

private:

    // the transformation parameters
    VectorF m_params;

    // disable copy constructor and assignment operator
    ImageTransform( const ImageTransform &x );
    ImageTransform &operator=( const ImageTransform &x );
};


} // end namespace sbl
#endif // _SBL_IMAGE_TRANSFORM_H_

