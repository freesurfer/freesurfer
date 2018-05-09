## 
## Makefile.am 
##

if HAVE_ITK_LIBS
ITK_THINGS=$(ITK_CFLAGS) -DHAVE_ITK_LIBS
else
ITK_THINGS=
endif

AM_CFLAGS=-I$(top_srcdir)/include -I$(top_srcdir)/include/dicom -O0
AM_CXXFLAGS=-I$(top_srcdir)/include -I$(top_srcdir)/include/dicom
AM_CPPFLAGS=-I$(top_srcdir)/include -I$(top_srcdir)/include/dicom -O0 $(ITK_THINGS) \
	-I/usr/pubsw/packages/itk/current/include/InsightToolkit/Utilities/vxl/core -I/usr/pubsw/packages/itk/current/include/InsightToolkit -I/usr/pubsw/packages/itk/current/include/InsightToolkit/Algorithms -I/usr/pubsw/packages/itk/current/include/InsightToolkit/BasicFilters -I/usr/pubsw/packages/itk/current/include/InsightToolkit/Common -I/usr/pubsw/packages/itk/current/include/InsightToolkit/IO -I/usr/pubsw/packages/itk/current/include/InsightToolkit/Numerics -I/usr/pubsw/packages/itk/current/include/InsightToolkit/Numerics/Statistics -I/usr/pubsw/packages/itk/current/include/InsightToolkit/Review -I/usr/pubsw/packages/itk/current/include/InsightToolkit/Review/Statistics -I/usr/pubsw/packages/itk/current/include/InsightToolkit/SpatialObject -I/usr/pubsw/packages/itk/current/include/InsightToolkit/Utilities/vxl/vcl
bin_PROGRAMS = dngtester cpptester

dngtester_SOURCES=dngtester.c
dngtester_LDADD= $(top_builddir)/fsgdf/libfsgdf.a \
  $(addprefix $(top_builddir)/, $(LIBS_MGH))
dngtester_LDFLAGS=$(OS_LDFLAGS)

cpptester_SOURCES=cpptester.cpp
cpptester_LDADD= $(top_builddir)/fsgdf/libfsgdf.a \
  $(addprefix $(top_builddir)/, $(LIBS_MGH))
cpptester_LDFLAGS=$(OS_LDFLAGS)

install:

# Our release target. Include files to be excluded here. They will be
# found and removed after 'make install' is run during the 'make
# release' target.
EXCLUDE_FILES=""
include $(top_srcdir)/Makefile.extra
