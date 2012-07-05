##
## Makefile.am 
##

foodir=$(prefix)/docs/xml
foo_DATA=mris_pmake.help.xml

bin_PROGRAMS = mris_pmake
mris_pmake_SOURCES = \
	abs2rel.c \
	asynch.cpp \
	c_label.cpp \
	C_mpmProg.cpp \
	c_SMessage.cpp \
	c_SSocket.cpp \
	c_surface.cpp \
	c_vertex.cpp \
	dijkstra.cpp \
	env.cpp \
	general.cpp \
	help.cpp \
	mris_pmake.cpp \
	rel2abs.c \
	scanopt.cpp \
	surface.cpp \
	C_mpmOverlay.cpp \
	C_mpmOverlay.h \
	legacy.cpp \
	legacy.h

if BUILDOPENCL # Include OpenCL version of algorithms
    AM_CFLAGS=-I$(top_srcdir)/include
    AM_CXXFLAGS=-I$(top_srcdir)/include \
	-I$(top_srcdir)/opencl_algorithms/ $(OPENCL_CFLAGS)
    mris_pmake_LDADD = $(addprefix $(top_builddir)/, $(LIBS_MGH)) \
	../opencl_algorithms/libopencl_algorithms.a $(OPENCL_LIBS)
    mris_pmake_LDFLAGS = $(OS_LDFLAGS)
else # No OpenCL present
    AM_CFLAGS=-I$(top_srcdir)/include
    AM_CXXFLAGS=-I$(top_srcdir)/include -I$(top_srcdir)/opencl_algorithms/ 
    mris_pmake_LDADD = $(addprefix $(top_builddir)/, $(LIBS_MGH)) \
	../opencl_algorithms/libopencl_algorithms.a
    mris_pmake_LDFLAGS = $(OS_LDFLAGS)
endif

TESTS=$(top_builddir)/scripts/help_xml_validate test_mris_pmake

EXTRA_DIST=asynch.h c_label.h C_mpmProg.h c_SMessage.h c_SSocket.h \
	c_surface.h c_vertex.h dijkstra.h env.h general.h help.h pathconvert.h \
	pstream.h scanopt.h surface.h $(foo_DATA) $(BUILT_SOURCES) test_mris_pmake

# Our release target. Include files to be excluded here. They will be
# found and removed after 'make install' is run during the 'make
# release' target.
EXCLUDE_FILES=
include $(top_srcdir)/Makefile.extra

BUILT_SOURCES=mris_pmake.help.xml.h

