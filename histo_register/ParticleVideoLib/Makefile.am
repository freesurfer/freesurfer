#
# Makefile.am
#

noinst_LIBRARIES=libpvl.a

libpvl_a_CXXFLAGS=-I$(top_srcdir)/histo_register/ParticleVideoLib/include \
	$(OPENCV_CXXFLAGS) -I$(top_srcdir)/histo_register/SimpleBaseLib/include
libpvl_a_CPPFLAGS=-DUSE_CDT -DUSE_OPENCV
libpvl_a_SOURCES=\
	src/Occlusion.cc \
	src/SimpleParticleBuild.cc \
	src/SimpleParticleData.cc \
	src/SimpleParticleOptimize.cc \
	src/SparseSystem.cc \
	src/VarMotion.cc \
	src/VarMotionMultiRes.cc \
	src/VarMotionUtil.cc \
	include/pvl/Occlusion.h \
	include/pvl/SimpleParticleBuild.h \
	include/pvl/SimpleParticleData.h \
	include/pvl/SimpleParticleOptimize.h \
	include/pvl/SparseSystem.h \
	include/pvl/VarMotionUtil.h \
	include/pvl/VarMotionMultiRes.h \
	include/pvl/VarMotion.h

