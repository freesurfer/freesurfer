#
# Makefile.am
#

AM_CXXFLAGS=-I$(top_srcdir)/include -I$(top_srcdir)/include/dicom
AM_LDFLAGS=-L$(top_srcdir)/lib

check_PROGRAMS=mnitest i2rtest icotest surftest extest mghxform inftest difftool
TESTS=mnitest i2rtest icotest extest mghxform inftest

mnitest_SOURCES=mnitest.cpp
i2rtest_SOURCES=i2rtest.cpp
icotest_SOURCES=icotest.cpp
surftest_SOURCES=surftest.cpp
extest_SOURCES=extest.cpp
mghxform_SOURCES=mghxform.cpp
inftest_SOURCES=inftest.cpp
difftool_SOURCES=difftool.cpp

LDADD=$(addprefix $(top_srcdir)/, $(LIBS_MGH))

