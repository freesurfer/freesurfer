#ifndef _PVL_VAR_MOTION_UTIL_H_
#define _PVL_VAR_MOTION_UTIL_H_
#include <sbl/core/Config.h>
#include <sbl/math/Vector.h>
#include <sbl/image/MotionField.h>
using namespace sbl;
namespace pvl {


/*! \file VarMotionUtil.h
    \brief The VarMotionUtil module provides various small utility functions used
    by the VarMotion and VarMotionMultiRes modules.
*/


//-------------------------------------------
// MISC UTILITY FUNCTIONS
//-------------------------------------------


/// returns vector of scales for multi-res optimization
VectorF computeScales( float scaleFactor, float minScale );


/// allocate and build mapping of flow components to system variables;
/// each pixel in mask corresponds to 2 system variables (du and dv);
/// returns count of pixels in mask
int buildIndexMaps( const ImageGrayU &mask, aptr<ImageGrayI> &uIndex, aptr<ImageGrayI> &vIndex );


/// create a motion field of the desired size;
/// if init is specified, use it (resize as appropriate)
aptr<MotionField> createInitMotion( int width, int height, const MotionField *init );


/// compute local smoothness factor (higher where image gradient is low)
aptr<ImageGrayF> localSmoothness( const ImageGrayF &src, float flatSigma, float flatSmoothness );


//-------------------------------------------
// DISCRETE DERIVATIVES
//-------------------------------------------


/// discrete derivatives
float dx( const ImageGrayF &img, int x, int y );
float dy( const ImageGrayF &img, int x, int y );


/// apply discrete derivative to image
aptr<ImageGrayF> dx( const ImageGrayF &img );
aptr<ImageGrayF> dy( const ImageGrayF &img );


/// apply discrete derivative to each channel
Array<ImageGrayF> dx( const Array<ImageGrayF> &chan );
Array<ImageGrayF> dy( const Array<ImageGrayF> &chan );


//-------------------------------------------
// HANDLE MULTI-CHANNEL IMAGES
//-------------------------------------------


/// displays each defined channel
void dispChannels( const Array<ImageGrayF> &chan );


/// scale and blur each channel
Array<ImageGrayF> shrinkChannels( const Array<ImageGrayF> &chan, int scaledWidth, int scaledHeight, float blurSigma );


/// given a video frame, extract a set of channels used for the data term
Array<ImageGrayF> extractChannels( const ImageColorU &img, Config &varConf );


/// apply a Gaussian blur to each channel
Array<ImageGrayF> blurChannels( const Array<ImageGrayF> &chan, float blurSigma );


/// allocate multi-channel image
Array<ImageGrayF> allocChannels( int width, int height, int count, float init );


//-------------------------------------------
// VARIATIONAL OBJECTIVE FUNCTION
//-------------------------------------------


/// this is the robust norm used in Brox et al. (denoted as Psi)
float broxDist( float diffSqd );


/// derivative of robust norm w.r.t. squared diff
float broxDistDeriv( float diffSqd );


/// the variational objective function 
float varEnergy( const MotionField &mf, const ImageGrayU &mask, const Array<ImageGrayF> &srcChanBlurred, const Array<ImageGrayF> &destChanBlurred, 
                 Config &varConf, bool verbose, float &dataEnergyRet, float &smoothnessEnergyRet,
                 ImageGrayF *dataEnergyMap, ImageGrayF *smoothnessEnergyMap );


/// the variational objective function
float varEnergy( const MotionField &mf, const ImageGrayU &mask, const ImageColorU &src, const ImageColorU &dest, 
                 Config &varConf, float &dataEnergy, float &smoothnessEnergy,
                 ImageGrayF *dataEnergyMap, ImageGrayF *smoothnessEnergyMap );


/// the variational objective function
float varEnergy( const MotionField &mf, const ImageGrayU &mask, 
                 const Array<ImageGrayF> &srcChanBlurred, const Array<ImageGrayF> &destChanBlurred,
                 Config &varConf, ImageGrayF &energyMap );


} // end namespace pvl
#endif // _PVL_VAR_MOTION_UTIL_H_

