#ifndef _SBL_IMAGE_UTIL_H_
#define _SBL_IMAGE_UTIL_H_
#include <sbl/core/Pointer.h>
#include <sbl/other/TaggedFile.h>
#include <sbl/image/Image.h>
namespace sbl {


/*! \file ImageUtil.h
    \brief The ImageUtil module provides various functions for working with image objects,
    including filtering, converting, and loading/saving.  Most of these functions are
    simple wrappers for OpenCV functions.  For spatial image transformations, see 
    the ImageTransform module.
*/


// register commands, etc. defined in this module
void initImageUtil();


//-------------------------------------------
// IMAGE CONVERSION
//-------------------------------------------


/// convert color image to gray image
aptr<ImageGrayU> toGray( const ImageColorU &img );


/// convert gray image to color image
aptr<ImageColorU> toColor( const ImageGrayU &img );


/// convert 8-bit image to float image
aptr<ImageGrayF> toFloat( const ImageGrayU &input, float scaleFactor );
aptr<ImageColorF> toFloat( const ImageColorU &input, float scaleFactor );


/// convert float image to 8-bit image, automatically scaling values 
aptr<ImageGrayU> toUChar( const ImageGrayF &img );


/// convert float image to 8-bit image, using a fixed scale factor
aptr<ImageGrayU> toUChar( const ImageGrayF &input, float scaleFactor );
aptr<ImageColorU> toUChar( const ImageColorF &input, float scaleFactor );


//-------------------------------------------
// IMAGE FILTERS
//-------------------------------------------


/// blur using box filter
template <typename ImageType> aptr<ImageType> blurBox( const ImageType &input, int boxSize );


/// blur using Gaussian filter
template <typename ImageType> aptr<ImageType> blurGauss( const ImageType &input, float sigma );


/// apply median filter (set each pixel to median of neighbors)
template <typename ImageType> aptr<ImageType> median( const ImageType &input, int apertureSize );


/// image x gradient
template <typename ImageType> aptr<ImageGrayF> xGrad( const ImageType &input, int apertureSize );


/// image y gradient
template <typename ImageType> aptr<ImageGrayF> yGrad( const ImageType &input, int apertureSize );


/// compute magnitude of each gradient vector
aptr<ImageGrayF> gradientMagnitude( const ImageGrayU &input, int apertureSize );


/// apply threshold to image values
template <typename ImageType> aptr<ImageType> threshold( const ImageType &input, float thresh, bool invert );


/// invert a floating point image, assuming values in [0, 1]
aptr<ImageGrayF> invert( const ImageGrayF &input );


/// apply a box filter and then threshold
aptr<ImageGrayU> blurBoxAndThreshold( const ImageGrayU &input, int boxSize, int thresh );
aptr<ImageGrayU> blurBoxAndThreshold( const ImageGrayF &input, int boxSize, float thresh );


/// multiply each pixel value by the specified factor
void multiply( const ImageGrayF &input, float factor, ImageGrayF &output );


/// perform local brightness normalization
aptr<ImageColorU> normalizeLocalBrightness( const ImageColorU &input, int windowSize );
aptr<ImageGrayU> normalizeLocalBrightness( const ImageGrayU &input, int windowSize );


//-------------------------------------------
// IMAGE STATISTICS
//-------------------------------------------


/// compute mean pixel value
float mean( const ImageGrayU &img );
float mean( const ImageGrayF &img );
float mean( const ImageColorU &img );


/// compute mean pixel value of each color channel
void channelMean( const ImageColorU &img, float &rMean, float &gMean, float &bMean );


/// compute mean absolute difference between pixel values
float meanAbsDiff( const ImageColorU &img1, const ImageColorU &img2 );
float meanAbsDiff( const ImageGrayU &img1, const ImageGrayU &img2, int xBorder, int yBorder );


/// compute mutual info between a pair of images
float mutualInfo( const ImageGrayU &img1, const ImageGrayU &img2, int xBorder, int yBorder, int bucketCount );


/// compute bounds of the non-zero mask region
void maskBounds( const ImageGrayU &mask, int &xMin, int &xMax, int &yMin, int &yMax );


/// count number of non-zero entries in mask
int maskCount( const ImageGrayU &mask );


/// compute min/mean/max pixel values
void imageStats( const ImageGrayU &img, int &min, float &mean, int &max );
void imageStats( const ImageGrayF &img, float &min, float &mean, float &max );


/// compute histogram of image pixel values
VectorI imageHistogram( const ImageGrayU &image, int xMin, int xMax, int yMin, int yMax );


//-------------------------------------------
// IMAGE FILE I/O
//-------------------------------------------


/// save image to file (format determined by extension)
template <typename ImageType> void saveImage( const ImageType &img, const String &fileName );


/// load image from file (format determined by extension)
template <typename ImageType> aptr<ImageType> load( const String &fileName );


//-------------------------------------------
// FLOOD FILL
//-------------------------------------------


/// fills starting at the given point, all values within the given range;
/// the fill value must *not* be within the fill range
int floodFill( ImageGrayU &img, int minRegionColor, int maxRegionColor, int fillValue, 
               int x, int y, float *xCentRet = NULL, float *yCentRet = NULL, 
               int *xMinRet = NULL, int *xMaxRet = NULL, 
               int *yMinRet = NULL, int *yMaxRet = NULL );


/// keep only mask components that have a pixel count within the given range (assumes mask has only 0 and 255 values)
void filterMaskComponents( ImageGrayU &mask, int minSize, int maxSize );


/// fill holes in mask smaller than given size (assumes mask has only 0 and 255 values)
void fillMaskHoles( ImageGrayU &mask, int maxSize );


//-------------------------------------------
// OTHER IMAGE UTILS
//-------------------------------------------


/// returns true if specified line intersects non-zero mask values;
/// (assumes points are inside image bounds)
bool lineIntersects( const ImageGrayU &mask, float x1, float y1, float x2, float y2 );


/// returns true if rectange intersects non-zero mask values;
/// (assumes points are inside image bounds)
bool rectIntersects( const ImageGrayU &mask, int x1, int y1, int x2, int y2 );


/// draw the an iso-contour of the mask on the given output image
void drawMaskBoundary( ImageColorU &output, const ImageGrayU &mask, int thresh, int r, int g, int b, int xBorder = 0, int yBoder = 0 );


/// join images horizontally
aptr<ImageColorU> joinHoriz( const ImageColorU &img1, const ImageColorU &img2 );
aptr<ImageGrayU> joinHoriz( const ImageGrayU &img1, const ImageGrayU &img2 );


/// blend two images, apply alpha weight to each image (alpha values typically in [0, 1])
aptr<ImageColorU> blend( const ImageColorU &image1, const ImageColorU &image2, float alpha1, float alpha2 );


/// show left half of image 1 with right half of image 2
aptr<ImageColorU> splitScreen( const ImageColorU &img1, const ImageColorU &img2 );


/// display an image (if using GUI)
void dispImage( const ImageColorU &img, const String &label = "" );
void dispImage( const ImageGrayU &img, const String &labe = "" );


} // end namespace sbl
#endif // _SBL_IMAGE_UTIL_H_

