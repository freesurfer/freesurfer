#
# Makefile.am
#

AM_CXXFLAGS=-I$(top_srcdir)/include -I$(top_srcdir)/include/dicom
AM_LDFLAGS=-L$(top_srcdir)/lib

check_PROGRAMS=mnitest i2rtest icotest surftest extest
TESTS=mnitest i2rtest icotest extest

mnitest_SOURCES=mnitest.cpp
i2rtest_SOURCES=i2rtest.cpp
icotest_SOURCES=icotest.cpp
surftest_SOURCES=surftest.cpp
extest_SOURCES=extest.cpp

LDADD=$(addprefix $(top_srcdir)/, $(LIBS_MGH))

