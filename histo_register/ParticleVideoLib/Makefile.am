#
# Makefile.am
#

noinst_LIBRARIES=libpvl.a

libpvl_a_CXXFLAGS=-Iinclude/ $(OPENCV_CXXFLAGS) -I$(top_srcdir)/histo_register/SimpleBaseLib/include
libpvl_a_CPPFLAGS=-DUSE_CDT
#libpvl_a_CXXFLAGS=$(AM_CXXFLAGS)
libpvl_a_SOURCES=\
	src/Occlusion.cc \
	src/SimpleParticleBuild.cc \
	src/SimpleParticleData.cc \
	src/SimpleParticleOptimize.cc \
	src/SparseSystem.cc \
	src/VarMotion.cc \
	src/VarMotionMultiRes.cc \
	src/VarMotionUtil.cc

