#ifndef _SBL_MOTION_FIELD_UTIL_H_
#define _SBL_MOTION_FIELD_UTIL_H_
#include <sbl/core/Pointer.h>
#include <sbl/math/Vector.h>
#include <sbl/image/MotionField.h>
namespace sbl {


/*! \file MotionFieldUtil.h
    \brief The MotionFieldUtil module includes functions for manipulating and analyzing
    MotionField objects (2D vector fields representing pixel correspondences between 
    a pair of images).
*/


//-------------------------------------------
// MOTION FIELD STATISTICS
//-------------------------------------------


/// get min/mean/max of u and v components 
void motionFieldStats( const MotionField &mf, float &uMin, float &vMin, float &uMean, float &vMean, float &uMax, float &vMax );


/// display statistics for motion field
void dispMotionStats( int indent, const MotionField &mf );


/// compute mean magnitude of motion vectors
float motionFieldMag( const MotionField &mf );


//-------------------------------------------
// MOTION FIELD VISUALIZATION
//-------------------------------------------

/// create color visualization of motion field; different colors represent different motions;
/// gray corresponds to no motion; black indicates occluded regions
aptr<ImageColorU> colorizeMotion( const MotionField &mf );


/// display motion field visualization
void dispMotion( const MotionField &mf, const String &caption = "" );


//-------------------------------------------
// MOTION FIELD FILE I/O
//-------------------------------------------


/// load motion field from binary data file
aptr<MotionField> loadMotionField( const String &fileName );


/// load motion field from floating-point image file
aptr<MotionField> loadMotionFieldImage( const String &fileName );


/// save motion field to binary data file
void saveMotionField( const MotionField &mf, const String &fileName );


/// save motion field to floating-point image file
void saveMotionFieldImage( const MotionField &mf, const String &fileName );


/// save motion field to text file
void saveMotionFieldText( const MotionField &mf, const String &fileName );


//-------------------------------------------
// MISC. MOTION FIELD UTILS
//-------------------------------------------


/// compute motion field divergence (du/dx + dv/dy)
aptr<ImageGrayF> motionDivergence( const MotionField &mf );


/// compute motion field gradient magnitude squared [(du/dx)^2 + (du/dy)^2 + (dv/dx)^2 + (dv/dy)^2]
aptr<ImageGrayF> motionGradMagSqd( const MotionField &mf );


/// compute mean motion field gradient magnitude
float meanGradMag( const MotionField &mf );


/// interpolate (triangulate) offsets to create motion field
aptr<MotionField> computeDenseMotion( const VectorF &x, const VectorF &y, const VectorF &u, const VectorF &v, 
                                      int width, int height );


/// compute inverse motion field using point sample triangulation;
/// sampleStep determines distance between samples
aptr<MotionField> invertMotionField( const MotionField &mf, const ImageGrayU &mask, int sampleStep );


} // end namespace sbl
#endif // _SBL_MOTION_FIELD_UTIL_H_

