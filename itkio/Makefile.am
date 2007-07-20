#
# Makefile.am for libitkio.a
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

noinst_LIBRARIES= libitkio.a

libitkio_a_SOURCES= \
  itkMGHImageIO.h itkMGHImageIO.cxx \
  itkMGHImageIOFactory.h itkMGHImageIOFactory.cxx 

am_libitkio_a_OBJECTS=$(libitkio_a_SOURCES:%.cxx=%.o)

endif
#end itk

EXTRA_DIST=$(wildcard *.h)

# Our release target. Include files to be excluded here. They will be
# found and removed after 'make install' is run during the 'make
# release' target.
EXCLUDE_FILES=
include $(top_srcdir)/Makefile.extra

	