#ifndef _VAR_CORRES_3D_H_
#define _VAR_CORRES_3D_H_
#include <sbl/core/Array.h>
#include <sbl/core/Config.h>
#include <sbl/core/Pointer.h>
#include "registration/CorresField3D.h"
using namespace sbl;
namespace hb {


/*! \file VarCorres3D.cc
    \brief The VarCorres3D module provides an algorithm for estimating a 3D mapping 
	(analogous to optical flow) from one image volume to another image volume.
*/


/// estimate a 3D mapping (analogous to optical flow) from one image volume to another image volume
void varCorres3D( const Array<ImageGrayU> &src, const Array<ImageGrayU> &dest, Array<CorresField3D> &cfSeq, int scaleFactor );


} // end namespace hb
#endif // _VAR_CORRES_3D_H_
