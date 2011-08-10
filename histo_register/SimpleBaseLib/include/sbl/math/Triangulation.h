#ifndef _SBL_TRIANGULATION_H_
#define _SBL_TRIANGULATION_H_
#ifdef USE_CDT
#include <sbl/math/Vector.h>
#include <sbl/math/Matrix.h>
#include <sbl/core/Pointer.h>
#include <sbl/image/Image.h>
namespace sbl {


/*! \file Triangulation.h
    \brief The Triangulation module includes functions for triangulating points
    using the CDT algorithm.
*/


// returns matrix of edges; each row is (x1, y1, x2, y2)
aptr<MatrixF> triangulate( const MatrixF &points );


// returns matrix of edges: each row is pair (i1, i2) of point indices
aptr<MatrixI> triangulateIndex( const MatrixF &points );


// use triangulation to interpolate values across an image
aptr<ImageGrayF> triangulationInterpolation( const MatrixF &points, VectorF &values, int width, int height );


} // end namespace sbl
#endif // USE_CDT
#endif // _SBL_TRIANGULATION_H_

