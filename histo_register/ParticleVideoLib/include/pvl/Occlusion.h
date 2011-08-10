#ifndef _PVL_OCCLUSION_H_
#define _PVL_OCCLUSION_H_
#include <sbl/core/Config.h>
#include <sbl/image/MotionField.h>
using namespace sbl;
namespace pvl {


/*! \file Occlusion.h
    \brief The Occlusion module provides code for detecting/labelling occlusions 
    in a MotionField object.
*/


/// run bilateral filter loop with occlusion detection
void filterMotionLoop( MotionField &mf, const ImageGrayF &src, const ImageGrayF &dest, Config &varConf, bool fast );


} // end namespace pvl
#endif // _PVL_OCCLUSION_H_

