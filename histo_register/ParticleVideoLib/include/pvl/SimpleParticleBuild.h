#ifndef _SIMPLE_PARTICLE_BUILD_H_
#define _SIMPLE_PARTICLE_BUILD_H_
#include <sbl/image/MotionField.h>
#include <pvl/SimpleParticleData.h>
namespace pvl {


/*! \file SimpleParticleBuild.h
    \brief The SimpleParticleBuild module gives the top-level functions for constructing
    a particle set and clustering a particle set.
*/


/// extend current particle set by one frame;
/// to build a particle set, call this function on consecutive frames (starting with frameIndex 0);
/// mf maps frameIndex to frameIndex + 1;
/// mfPrev maps frameIndex - 1 to frameIndex (or NULL if frameIndex == 0)
void simpleParticleFrame( SimpleParticleSet &particleSet, int frameIndex, const ImageColorU &frame, 
                          const MotionField *mfPrev, const MotionField &mf, Config &pfConf );


/// create particle clusters using a particle-pair distance measure
void buildParticleClusters( SimpleParticleSet &particleSet, int width, int height, Config &spConf, bool visualize );


} // end namespace pvl
#endif // _SIMPLE_PARTICLE_BUILD_H_

