#ifndef _SBL_IMAGE_SEQ_UTIL_H_
#define _SBL_IMAGE_SEQ_UTIL_H_
#include <sbl/core/Array.h>
#include <sbl/image/Image.h>
namespace sbl {


/*! \file ImageSeqUtil.h
    \brief The ImageSeqUtil module includes functions for manipulating image
    sequences.  The functions assume that all the images in a sequence have 
    the same dimensions.
*/


// commonly used image sequence types
typedef Array<ImageGrayU> ImageGrayUSeq;
typedef Array<ImageGrayF> ImageGrayFSeq;
typedef Array<ImageGrayI> ImageGrayISeq;


/// create a sequence of images
void initImageSeq( ImageGrayUSeq &seq, int width, int height, int length, bool clear, int clearValue );
void initImageSeq( ImageGrayFSeq &seq, int width, int height, int length, bool clear, float clearValue );
void initImageSeq( ImageGrayISeq &seq, int width, int height, int length, bool clear, int clearValue );


/// resize (in x, y) a sequence of images
void resizeSeqXY( ImageGrayUSeq &seq, int newWidth, int newHeight );
void resizeSeqXY( ImageGrayFSeq &seq, int newWidth, int newHeight );


/// resize (in z) a sequence of images
void resizeSeqZ( ImageGrayUSeq &inSeq, ImageGrayUSeq &outSeq, int newLength );
void resizeSeqZ( ImageGrayFSeq &inSeq, ImageGrayFSeq &outSeq, int newLength );


/// blur (in x, y) a sequence of images
void blurGaussSeqXY( ImageGrayUSeq &seq, float sigma );
void blurGaussSeqXY( ImageGrayFSeq &seq, float sigma );


/// blur (in z) a sequences of images
void blurGaussSeqZ( const ImageGrayUSeq &inSeq, ImageGrayUSeq &outSeq, float sigma );
void blurGaussSeqZ( const ImageGrayFSeq &inSeq, ImageGrayFSeq &outSeq, float sigma );


/// perform bilinear interpolation to find value in image sequence
float interp( const ImageGrayUSeq &seq, float x, float y, float z );
float interp( const ImageGrayFSeq &seq, float x, float y, float z );


/// transpose the Y and Z axes in the image sequence
void transposeYZ( const ImageGrayUSeq &inSeq, ImageGrayUSeq &outSeq );


/// fills starting at the given point, all values within the given range;
/// the fill value must *not* be within the fill range
int floodFillXYZ( const ImageGrayUSeq &seq, int minRegionColor, int maxRegionColor, int fillValue, 
                  int x, int y, int z,
                  float *xCentRet = NULL, float *yCentRet = NULL, float *zCentRet = NULL, 
                  int *xMinRet = NULL, int *xMaxRet = NULL, 
                  int *yMinRet = NULL, int *yMaxRet = NULL,
                  int *zMinRet = NULL, int *zMaxRet = NULL );


} // end namespace sbl
#endif // _SBL_IMAGE_SEQ_UTIL_H_

