#ifndef _PVL_SIMPLE_PARTICLE_OPTIMIZE_H_
#define _PVL_SIMPLE_PARTICLE_OPTIMIZE_H_
#include <pvl/SimpleParticleData.h>
#include <sbl/math/Matrix.h>
using namespace sbl;
namespace pvl {


/*! \file SimpleParticleOptimize.h
    \brief The SimpleParticleOptimize module provides functions for optimizing the
    locations of a set of particles.
*/


/// optimizes particle positions using blurred channel values
void simpleParticleOptimize( SimpleParticleSet &particleSet, int frameIndex, 
                             const Array<ImageGrayF> &chan, Config &spConf );


/// returns a matrix in which each row is a (particleIndex1, particleIndex2) indicating a link between a pair of particle;
/// if activePrev set, requires that particle be active in previous frame
aptr<MatrixI> createLinks( const SimpleParticleSet &particleSet, int frameIndex, bool activePrev );


} // end namespace pvl
#endif // _PVL_SIMPLE_PARTICLE_OPTIMIZE_H_

