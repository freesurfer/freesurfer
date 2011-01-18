#
# Makefile.am for libopencl_algorithms.a
#

noinst_LIBRARIES=libopencl_algorithms.a

if BUILDOPENCL
# If OpenCL was detected on the system, include all of the OpenCL code.
# Otherwise, just build versions of the code that do not require
# OpenCL (reference algorithms).
libopencl_algorithms_a_SOURCES = oclDijkstraKernel.cpp oclDijkstraKernel.h \
	oclCommon.cpp oclCommon.h
libopencl_algorithms_a_CXXFLAGS = $(AM_CXXFLAGS) $(OPENCL_CFLAGS)
else
libopencl_algorithms_a_SOURCES = oclDijkstraKernel.cpp
endif

EXTRA_DIST = dijkstra.cl dijkstra.cl.h make_cl_header.sh

# Our release target. Include files to be excluded here. They will be
# found and removed after 'make install' is run during the 'make
# release' target.
EXCLUDE_FILES=
include $(top_srcdir)/Makefile.extra
