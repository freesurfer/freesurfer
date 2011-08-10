#ifndef _PVL_VAR_MOTION_MULTI_RES_H_
#define _PVL_VAR_MOTION_MULTI_RES_H_
#include <sbl/core/Array.h>
#include <sbl/core/Config.h>
#include <sbl/image/MotionField.h>
using namespace sbl;
namespace pvl {


/*! \file VarMotionMultiRes.h
    \brief The VarMotionMultiRes module provides a variation of the varMotion algorithm
    in which the flow field is represented at a lower resolution than the image data
    (allowing faster optical flow estimation).
*/


/// computes an update (du, dv) to the flow field using a different resolution for the image vs. the motion field
void varMultiResIteration( MotionField &mf, const Array<ImageGrayF> &srcChanScaled, const Array<ImageGrayF> &destChanScaled, 
                           const ImageGrayU &maskScaled, Config &varConf );


} // end namespace pvl
#endif // _PVL_VAR_MOTION_MULTI_RES_H_

