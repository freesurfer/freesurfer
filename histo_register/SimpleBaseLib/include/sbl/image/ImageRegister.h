#ifndef _SBL_IMAGE_REGISTER_H_
#define _SBL_IMAGE_REGISTER_H_
#include <sbl/image/ImageTransform.h>
namespace sbl {


/*! \file ImageRegister.h
    \brief The ImageRegister module provides functions for estimating a transformation
    that aligns one image to another.
*/


/// computes the mean-abs image difference given an image transformation;
/// a step value greater than one allows faster (less accurate) optimization by ignoring some pixels;
/// the xBorder and yBorder specify areas to be ignored for the objective function;
/// the offsetBound value specifies the largest allowed translation between the images
double evalImageTransform( const ImageTransform &transform, const ImageGrayU &src, const ImageGrayU &dest, int step, int border, const ImageGrayU *srcMask = NULL, const ImageGrayU *destMask = NULL, bool interp = false, bool verbose = false );


/// registers a pair of images using to minimize mean-abs difference;
/// a step value greater than one allows faster (less accurate) optimization by ignoring some pixels;
/// the xBorder and yBorder specify areas to be ignored for the objective function;
/// the offsetBound value specifies the largest allowed translation between the images
aptr<ImageTransform> registerUsingImageTransform( const ImageGrayU &src, const ImageGrayU &dest, int transformParamCount, int step, int xBorder, int yBorder, float offsetBound, const ImageTransform *initTransform = NULL, const ImageGrayU *srcMask = NULL, const ImageGrayU *destMask = NULL, bool interp = false );


} // end namespace sbl
#endif // _SBL_IMAGE_REGISTER_H_

