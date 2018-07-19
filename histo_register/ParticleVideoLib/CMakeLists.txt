project(ParticleVideoLib)

set(SOURCES
  src/Occlusion.cc
  src/SimpleParticleBuild.cc
  src/SimpleParticleData.cc
  src/SimpleParticleOptimize.cc
  src/SparseSystem.cc
  src/VarMotion.cc
  src/VarMotionMultiRes.cc
  src/VarMotionUtil.cc
)

add_library(pvl STATIC ${SOURCES})