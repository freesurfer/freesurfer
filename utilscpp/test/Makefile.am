#
# Makefile.am
#

# we pick linker to be g++, since it knows where to find 
# the right libstdc++ and other
# files needed to resolve symbols in the library of C++ code
CCLD=g++

AM_CFLAGS=-I$(top_srcdir)/include \
	-I$(top_srcdir)/include/dicom

check_PROGRAMS=ctest

ctest_SOURCES=ctest.c
ctest_LDFLAGS=$(OS_LDFLAGS)

#segfaults: TESTS=ctest

LDADD=$(addprefix $(top_builddir)/, $(LIBS_MGH) utilscpp/libutilscpp.a)

# Our release target. Include files to be excluded here. They will be
# found and removed after 'make install' is run during the 'make
# release' target.
EXCLUDE_FILES=
include $(top_srcdir)/Makefile.extra
