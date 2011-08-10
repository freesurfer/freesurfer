#ifndef _VAR_CORRES_3D_UTIL_H_
#define _VAR_CORRES_3D_UTIL_H_
#include <sbl/core/Config.h>
#include <sbl/math/Vector.h>
#include "registration/CorresField3D.h"
#include "registration/ImageSetSeq.h"
using namespace sbl;
namespace hb {


/*! \file VarCorres3DUtil.cc
    \brief The VarCorres3DUtil module contains utility functions used by the VarCorres3D algorithm.
*/


//-------------------------------------------
// MISC UTILITY FUNCTIONS
//-------------------------------------------


/// returns vector of scales for multi-res optimization
VectorF varScaleSequence( float scaleFactor, float minScale );


/// build mapping of flow components to system variables;
/// each pixel in mask corresponds to 3 system variables (du, dv, and dw);
/// returns count of pixels in mask
int buildIndexMaps( int startIndex, const ImageGrayU &mask, ImageGrayI &uIndex, ImageGrayI &vIndex, ImageGrayI &wIndex );


/// returns config object with parameters for variational correspondence algorithm
aptr<Config> _varCorres3DConfig();


//-------------------------------------------
// DISCRETE DERIVATIVES
//-------------------------------------------


/// discrete derivatives
float dx( const ImageGrayF &img, int x, int y );
float dy( const ImageGrayF &img, int x, int y );
float dz( const Array<ImageGrayF> &imgSeq, int x, int y, int z );


/// apply discrete derivative to image
aptr<ImageGrayF> dx( const ImageGrayF &img );
aptr<ImageGrayF> dy( const ImageGrayF &img );
aptr<ImageGrayF> dz( const Array<ImageGrayF> &seq, int i );


/// apply discrete derivative to each channel
Array<ImageGrayF> dx( const Array<ImageGrayF> &chan );
Array<ImageGrayF> dy( const Array<ImageGrayF> &chan );
Array<ImageGrayF> dz( const ImageSetSeq &seq, int i );


/// compute discrete derivative of a single slice of a correspondence volume
aptr<ImageGrayF> dz( const Array<CorresField3D> &cfSeq, int i, int axis );


//-------------------------------------------
// VARIATIONAL OBJECTIVE FUNCTION
//-------------------------------------------


// this is the robust norm used in Brox et al. (denoted as Psi)
float psi( float diffSqd, bool robust );


// derivative of robust norm w.r.t. squared diff
float psiDeriv( float diffSqd, bool robust );


} // end namespace hb
#endif // _VAR_CORRES_3D_UTIL_H_
