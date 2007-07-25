#
# Makefile.am for libitkutils.a
#

#itk is required for building poistats
# add ITK
if HAVE_ITK_LIBS

ITK_THINGS=$(ITK_CFLAGS) -DHAVE_ITK_LIBS

AM_CPPFLAGS =	 -Wno-deprecated \
  $(X_CFLAGS) \
  $(VXL_CFLAGS) \
  -I$(top_srcdir)/include \
  $(ITK_THINGS) \
  -DHAVE_ZLIB

AM_LDFLAGS=

noinst_LIBRARIES= libitkutils.a

BSPLINES_SRCS = \
	bspline/itkBSplineKernelFunction.txx \
	bspline/itkBSplineKernelFunction.h \
	bspline/itkBSplineDataPointSetToImageFilter.txx \
  bspline/itkBSplineDataPointSetToImageFilter.h \
	bspline/itkPointSetToImageFilter.txx \
	bspline/itkPointSetToImageFilter.h

UI_SRCS = \
  ui/CommandParser.cxx \
  ui/FreeSurferExecutable.cxx

libitkutils_a_SOURCES= \
  $(BSPLINES_SRCS) \
  $(UI_SRCS) \
  itkOrientedVectorImage.h itkOrientedVectorImage.txx

# I don't think that I need this
#am_libitkutils_a_OBJECTS=$(libitkutils_a_SOURCES:%.cxx=%.o)

endif
#end itk

EXTRA_DIST=$(wildcard *.h) ui

# Our release target. Include files to be excluded here. They will be
# found and removed after 'make install' is run during the 'make
# release' target.
EXCLUDE_FILES=ui/CVS
include $(top_srcdir)/Makefile.extra
