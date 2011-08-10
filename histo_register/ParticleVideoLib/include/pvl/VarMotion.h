#ifndef _PVL_VAR_MOTION_H_
#define _PVL_VAR_MOTION_H_
#include <sbl/core/Array.h>
#include <sbl/core/Config.h>
#include <sbl/core/Pointer.h>
#include <sbl/image/MotionField.h>
using namespace sbl;
namespace pvl {


/*! \file VarMotion.h
    \brief The VarMotion module gives an implementation of a variational optical flow 
    algorithm similar to the Brox et al. algorithm.
*/


/// top-level variational (modified Brox et al.) motion estimation algorithm;
/// mask (if any) specifies where the algorithm should estimate flow vectors;
/// init (if any) provides an initial motion estimation;
/// scaleFactor (if not 1) is used to estimate flow at a lower resolution (for example, at half resolution if scaleFactor == 2)
aptr<MotionField> varMotion( const ImageColorU &src, const ImageColorU &dest, const String &configFileName, const MotionField *init, const ImageGrayU *mask = NULL, int scaleFactor = 1 );
aptr<MotionField> varMotion( const ImageGrayU &src, const ImageGrayU &dest, const String &configFileName, const MotionField *init, const ImageGrayU *mask = NULL, int scaleFactor = 1 );


/// computes an update (du, dv) to the flow field
void varIteration( MotionField &mf, const Array<ImageGrayF> &srcChan, const Array<ImageGrayF> &destChan, 
                   const ImageGrayU &mask, Config &varConf, bool fast );


} // end namespace pvl
#endif // _PVL_VAR_MOTION_H_

