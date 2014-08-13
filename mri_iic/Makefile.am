##
## Makefile.am 
##

AM_CFLAGS=-I$(top_srcdir)/include 
AM_CXXFLAGS=-I$(top_srcdir)/include $(VXL_CFLAGS)

check_PROGRAMS = testTCS \
	testMatrix \
	testMultiplication

testTCS_SOURCES=testTCS.cpp tensorCubicSmoothing.cpp matrix3d.cpp
testTCS_LDADD= $(addprefix $(top_builddir)/, $(LIBS_MGH))
testTCS_LDFLAGS=$(OS_LDFLAGS)

testMatrix_SOURCES=testMatrix.cpp matrix3d.cpp
testMatrix_LDADD= $(addprefix $(top_builddir)/, $(LIBS_MGH))
testMatrix_LDFLAGS=$(OS_LDFLAGS)

testMultiplication_SOURCES=testMultiplication.cpp matrix3d.cpp
testMultiplication_LDADD= $(addprefix $(top_builddir)/, $(LIBS_MGH))
testMultiplication_LDFLAGS=$(OS_LDFLAGS)

TESTS=testMultiplication testTCS

# Our release target. Include files to be excluded here. They will be
# found and removed after 'make install' is run during the 'make
# release' target.
EXCLUDE_FILES=
include $(top_srcdir)/Makefile.extra
