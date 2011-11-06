// Licensed under MIT license; see license.txt.

#include <sbl/image/ImageTransform.h>
#include <sbl/math/MathUtil.h>
#ifdef USE_OPENCV
    #include <opencv/cv.h>
    #include <opencv/highgui.h>
#endif
namespace sbl {


//-------------------------------------------
// IMAGE TRANSFORMATION
//-------------------------------------------


/// extract sub-image
// fix(clean): should unify with grayscale version
template<> aptr<ImageColorU> crop( const ImageColorU &input, int xMin, int xMax, int yMin, int yMax ) {
    assertDebug( xMin > 0 && xMax < input.width() );
    assertDebug( yMin > 0 && yMax < input.height() );
    int newWidth = xMax - xMin + 1;
    int newHeight = yMax - yMin + 1;
    aptr<ImageColorU> output( new ImageColorU( newWidth, newHeight ) );
    for (int y = yMin; y <= yMax; y++) 
        for (int x = xMin; x <= xMax; x++)
            for (int c = 0; c < 3; c++)
                output->data( x - xMin, y - yMin, c ) = input.data( x, y, c );
    return output;
}


/// extract sub-image
// fix(faster): use memcpy
template <typename ImageType> aptr<ImageType> crop( const ImageType &input, int xMin, int xMax, int yMin, int yMax ) {
    assertDebug( xMin > 0 && xMax < input.width() );
    assertDebug( yMin > 0 && yMax < input.height() );
    int newWidth = xMax - xMin + 1;
    int newHeight = yMax - yMin + 1;
    aptr<ImageType> output( new ImageType( newWidth, newHeight ) );
    for (int y = yMin; y <= yMax; y++) 
        for (int x = xMin; x <= xMax; x++)
            output->data( x - xMin, y - yMin ) = input.data( x, y );
    return output;
}
template aptr<ImageGrayU> crop( const ImageGrayU &input, int xMin, int xMax, int yMin, int yMax );
template aptr<ImageGrayF> crop( const ImageGrayF &input, int xMin, int xMax, int yMin, int yMax );


/// shrink or zoom image
template <typename ImageType> aptr<ImageType> resize( const ImageType &input, int newWidth, int newHeight, bool filter ) {
    aptr<ImageType> output( new ImageType( newWidth, newHeight ) );
#ifdef USE_OPENCV
    int width = input.width();
    cvResize( input.iplImage(), output->iplImage(), filter ? (newWidth > width ? CV_INTER_LINEAR : CV_INTER_AREA) : CV_INTER_NN );
#else
    fatalError( "not implemented" );
#endif
    return output;
}
template aptr<ImageGrayU>  resize( const ImageGrayU  &input, int newWidth, int newHeight, bool filter );
template aptr<ImageGrayF>  resize( const ImageGrayF  &input, int newWidth, int newHeight, bool filter );
template aptr<ImageColorU> resize( const ImageColorU &input, int newWidth, int newHeight, bool filter );


/// translate and scale an image
template <typename ImageType> aptr<ImageType> shiftScale( const ImageType &input, float xOffset, float yOffset, float xScale, float yScale, int outputWidth, int outputHeight ) {
    aptr<ImageType> output( new ImageType( outputWidth, outputHeight ) );
#ifdef USE_OPENCV
    CvMat *map = cvCreateMat( 2, 3, CV_32FC1 );
    cvmSet( map, 0, 0, xScale );
    cvmSet( map, 0, 1, 0);
    cvmSet( map, 0, 2, xOffset );
    cvmSet( map, 1, 0, 0);
    cvmSet( map, 1, 1, yScale);
    cvmSet( map, 1, 2, yOffset );
    cvWarpAffine( input.iplImage(), output->iplImage(), map, CV_INTER_CUBIC + CV_WARP_FILL_OUTLIERS, cvScalarAll( 255 ) );
    // fix(later): use cvGetQuadrangleSubPix?
    cvReleaseMat( &map );
#else
    fatalError( "not implemented" );
#endif
    return output;
}
template aptr<ImageGrayU>  shiftScale( const ImageGrayU  &input, float xOffset, float yOffset, float xScale, float yScale, int outputWidth, int outputHeight );
template aptr<ImageGrayF>  shiftScale( const ImageGrayF  &input, float xOffset, float yOffset, float xScale, float yScale, int outputWidth, int outputHeight );
template aptr<ImageColorU> shiftScale( const ImageColorU &input, float xOffset, float yOffset, float xScale, float yScale, int outputWidth, int outputHeight );


/// apply linear transformation to image
template <typename ImageType> aptr<ImageType> warpAffine( const ImageType &input, float xOffset, float yOffset, float x1, float y1, float x2, float y2, int outputWidth, int outputHeight, int fillColor ) {
    aptr<ImageType> output( new ImageType( outputWidth, outputHeight ) );
#ifdef USE_OPENCV
    CvMat *map = cvCreateMat( 2, 3, CV_32FC1 );
    cvmSet( map, 0, 0, x1 );
    cvmSet( map, 0, 1, x2 );
    cvmSet( map, 0, 2, xOffset );
    cvmSet( map, 1, 0, y1 );
    cvmSet( map, 1, 1, y2 );
    cvmSet( map, 1, 2, yOffset );
    cvWarpAffine( input.iplImage(), output->iplImage(), map, CV_INTER_CUBIC + CV_WARP_FILL_OUTLIERS, cvScalarAll( fillColor ) );
    // fix(later): use cvGetQuadrangleSubPix?
    cvReleaseMat( &map );
#else
    fatalError( "not implemented" );
#endif
    return output;
}
template aptr<ImageGrayU>  warpAffine( const ImageGrayU  &input, float xOffset, float yOffset, float x1, float y1, float x2, float y2, int outputWidth, int outputHeight, int fillColor );
template aptr<ImageGrayF>  warpAffine( const ImageGrayF  &input, float xOffset, float yOffset, float x1, float y1, float x2, float y2, int outputWidth, int outputHeight, int fillColor );
template aptr<ImageColorU> warpAffine( const ImageColorU &input, float xOffset, float yOffset, float x1, float y1, float x2, float y2, int outputWidth, int outputHeight, int fillColor );


/// flip image vertically (about horizontal axis)
template <typename ImageType> aptr<ImageType> flipVert( const ImageType &input ) {
    aptr<ImageType> output( new ImageType( input.width(), input.height() ) );
#ifdef USE_OPENCV
    cvConvertImage( input.iplImage(), output->iplImage(), CV_CVTIMG_FLIP );
#else
    fatalError( "not implemented" );
#endif
    return output;
}
template aptr<ImageGrayU>  flipVert( const ImageGrayU  &input );
template aptr<ImageGrayF>  flipVert( const ImageGrayF  &input );
template aptr<ImageColorU> flipVert( const ImageColorU &input );

/// flip image horizontally (about vertical axis)
template <typename ImageType> aptr<ImageType> flipHoriz( const ImageType &input ) {
    aptr<ImageType> output( new ImageType( input.width(), input.height() ) );
#ifdef USE_OPENCV
    cvFlip(input.iplImage(), output->iplImage(), 1);
#else
    fatalError( "not implemented" );
#endif
    return output;
}
template aptr<ImageGrayU>  flipHoriz( const ImageGrayU  &input );
template aptr<ImageGrayF>  flipHoriz( const ImageGrayF  &input );
template aptr<ImageColorU> flipHoriz( const ImageColorU &input );

// /// flip image horizontally (about vertical axis)
// // fix(clean): use opencv and templates
// aptr<ImageGrayU> flipHoriz( const ImageGrayU &input ) {
//     aptr<ImageGrayU> output( new ImageGrayU( input.width(), input.height() ) );
//     int width = input.width(), height = input.height();
//     for (int y = 0; y < height; y++) {
//         for (int x = 0; x < width; x++) {
//             output->data( x, y ) = input.data( width - x - 1, y );
//         }
//     }
//     return output;
// }


/// rotate image 180 degrees
// fix(clean): use opencv
template <typename ImageType> aptr<ImageType> rotate180( const ImageType &input ) {
    aptr<ImageType> output( new ImageType( input.width(), input.height() ) );
    int width = input.width(), height = input.height(), cc = input.channelCount();
    for (int y = 0; y < height; y++) {
        for (int x = 0; x < width; x++) {
            for (int c = 0; c< cc ; c++) {
                output->data( x, y, c ) = input.data( width - x - 1, height - y - 1, c );
            }
        }
    }
    return output;
}
template aptr<ImageGrayU>  rotate180( const ImageGrayU  &input );
template aptr<ImageGrayF>  rotate180( const ImageGrayF  &input );
template aptr<ImageColorU> rotate180( const ImageColorU &input );


/// rotate 90 degrees (counter-clockwise)
// fix(clean): use opencv
template <typename ImageType> aptr<ImageType> rotate90( const ImageType &input ) {
    int newWidth = input.height();
    int newHeight = input.width();
    int cc = input.channelCount();
    aptr<ImageType> output( new ImageType( newWidth, newHeight ) );
    for (int y = 0; y < newHeight; y++) {
        for (int x = 0; x < newWidth; x++) {
            for (int c = 0; c< cc ; c++) {
                output->data( x, y, c ) = input.data( y, newWidth - x , c );
            }
        }
    }
    return output;
}
template aptr<ImageGrayU>  rotate90( const ImageGrayU  &input );
template aptr<ImageGrayF>  rotate90( const ImageGrayF  &input );
template aptr<ImageColorU> rotate90( const ImageColorU &input );


/// rotate 270 degrees (counter-clockwise)
// fix(clean): use opencv
template <typename ImageType> aptr<ImageType> rotate270( const ImageType &input ) {
    int newWidth = input.height();
    int newHeight = input.width();
    int cc = input.channelCount();
    aptr<ImageType> output( new ImageType( newWidth, newHeight ) );
    for (int y = 0; y < newHeight; y++) {
        for (int x = 0; x < newWidth; x++) {
            for (int c = 0; c< cc ; c++) {
                output->data( x, y, c ) = input.data( newHeight - y, x , c );
            }
        }
    }
    return output;
}
template aptr<ImageGrayU>  rotate270( const ImageGrayU  &input );
template aptr<ImageGrayF>  rotate270( const ImageGrayF  &input );
template aptr<ImageColorU> rotate270( const ImageColorU &input );


/// rotate arbitrary amount (counter-clockwise)
template <typename ImageType> aptr<ImageType> rotate( const ImageType &input, float angleDegrees, int fillColor ) {
    int width = input.width(), height = input.height();
    float theta = angleDegrees * 3.14159f / 180.0f;
    float c = cosf( theta );
    float s = sinf( theta );
    float x1 = c;
    float y1 = s;
    float x2 = -s;
    float y2 = c;
    float xCenter = (float) width * 0.5f;
    float yCenter = (float) height * 0.5f;
    float xOffset = (1.0f - c) * xCenter + s * yCenter;
    float yOffset = -s * xCenter + (1.0f - c) * yCenter;
    int outputWidth = width;
    int outputHeight = height;
    if (fAbs( s ) > fAbs( c )) {
        outputWidth = height;
        outputHeight = width;
    }
    aptr<ImageType> output = warpAffine( input, xOffset, yOffset, x1, y1, x2, y2, outputWidth, outputHeight, fillColor );
    return output;
}
template aptr<ImageGrayU>  rotate( const ImageGrayU  &input, float angleDegrees, int fillColor );
template aptr<ImageGrayF>  rotate( const ImageGrayF  &input, float angleDegrees, int fillColor );
template aptr<ImageColorU> rotate( const ImageColorU &input, float angleDegrees, int fillColor );


//-------------------------------------------
// IMAGE TRANSFORM CLASS
//-------------------------------------------


/// create a translation transformation
ImageTransform::ImageTransform( float xOffset, float yOffset ) { 
    m_params.setLength( 2 ); 
    m_params[ 0 ] = xOffset; 
    m_params[ 1 ] = yOffset; 
}


/// create an affine transformation
ImageTransform::ImageTransform( float xOffset, float yOffset, float xScale, float yScale ) { 
    m_params.setLength( 6 ); 
    m_params[ 0 ] = xOffset; 
    m_params[ 1 ] = yOffset; 
    m_params[ 2 ] = xScale; 
    m_params[ 3 ] = 0; 
    m_params[ 4 ] = 0; 
    m_params[ 5 ] = yScale; 
}


/// load transformation parameters from file
ImageTransform::ImageTransform( File &file ) {
    m_params = file.readVector<float>();
}


/// save transformation parameters to file
void ImageTransform::save( File &file ) const {
    file.writeVector( m_params );
}


/// compute the inverse of the transformation
aptr<ImageTransform> ImageTransform::inverse() const {
    aptr<ImageTransform> invTransform;
    if (m_params.length() == 2) {
        VectorF invParams( 2 );
        invParams[ 0 ] = -m_params[ 0 ];
        invParams[ 1 ] = -m_params[ 1 ];
        invTransform.reset( new ImageTransform( invParams ) );
    } else if (m_params.length() == 6) {
        VectorF invParams( 6 );
        float a = m_params[ 2 ];
        float b = m_params[ 4 ];
        float c = m_params[ 3 ];
        float d = m_params[ 5 ];
        float factor = 1.0f / (a * d - b * c);
        invParams[ 2 ] = d * factor;
        invParams[ 4 ] = -b * factor;
        invParams[ 3 ] = -c * factor;
        invParams[ 5 ] = a * factor;
        invParams[ 0 ] = -(invParams[ 2 ] * m_params[ 0 ] + invParams[ 4 ] * m_params[ 1 ]);
        invParams[ 1 ] = -(invParams[ 3 ] * m_params[ 0 ] + invParams[ 5 ] * m_params[ 1 ]);
        invTransform.reset( new ImageTransform( invParams ) );
    } else {
        fatalError( "ImageTransform::inverse: not implemented" );
    }
    return invTransform;
}


/// map an image forward according to the transformation
aptr<ImageGrayU> ImageTransform::mapForward( const ImageGrayU &img, int outputWidth, int outputHeight, int fillColor ) const {
    if (m_params.length() == 6)
        return warpAffine( img, m_params[ 0 ], m_params[ 1 ], m_params[ 2 ], m_params[ 3 ], m_params[ 4 ], m_params[ 5 ], outputWidth, outputHeight, fillColor );    
    assertAlways( m_params.length() == 2)
    return shiftScale( img, m_params[ 0 ], m_params[ 1 ], 1, 1, outputWidth, outputHeight );
}


/// map an image forward according to the transformation
aptr<ImageColorU> ImageTransform::mapForward( const ImageColorU &img, int outputWidth, int outputHeight, int fillColor ) const {
    if (m_params.length() == 6)
        return warpAffine( img, m_params[ 0 ], m_params[ 1 ], m_params[ 2 ], m_params[ 3 ], m_params[ 4 ], m_params[ 5 ], outputWidth, outputHeight, fillColor );    
    assertAlways( m_params.length() == 2)
    return shiftScale( img, m_params[ 0 ], m_params[ 1 ], 1, 1, outputWidth, outputHeight );
}


/// map an image back using the inverse transformation
aptr<ImageGrayU> ImageTransform::mapBackward( const ImageGrayU &img, int outputWidth, int outputHeight, int fillColor ) const {
    aptr<ImageTransform> inverseTransform = inverse();
    return inverseTransform->mapForward( img, outputWidth, outputHeight, fillColor );
}


/// map a point forward according to the transformation
Point2 ImageTransform::mapForward( const Point2 &pt ) const {
    assertAlways( m_params.length() == 2 ); // fix(soon): handle affine xforms
    return Point2( pt.x + m_params[ 0 ], pt.y + m_params[ 1 ] );
}


/// map a point backward according to the inverse transformation
Point2 ImageTransform::mapBackward( const Point2 &pt ) const {
    assertAlways( m_params.length() == 2 ); // fix(soon): handle affine xforms
    return Point2( pt.x - m_params[ 0 ], pt.y - m_params[ 1 ] );
}


/// display transformation parameters
void ImageTransform::display( int indent ) {
    if (m_params.length() == 2) {
        disp( indent, "x: %f, y: %f", m_params[ 0 ], m_params[ 1 ] );
    } else if (m_params.length() == 6) {
        disp( indent, "x: %f, y: %f, x1: %f, y1: %f, x2: %ff, y2: %f", m_params[ 0 ], m_params[ 1 ], m_params[ 2 ], m_params[ 3 ], m_params[ 4 ], m_params[ 5 ] );
    } else {
        fatalError( "ImageTransform::display: not implemented" );
    } 
}


} // end namespace sbl

