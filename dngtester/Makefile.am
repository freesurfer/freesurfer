## 
## Makefile.am 
##

if HAVE_ITK_LIBS
ITK_THINGS=$(ITK_CFLAGS) -DHAVE_ITK_LIBS
else
ITK_THINGS=
endif

AM_CFLAGS=-I$(top_srcdir)/include -I$(top_srcdir)/include/dicom -O0
AM_CXXFLAGS=-I$(top_srcdir)/include -I$(top_srcdir)/include/dicom -O0 $(ITK_THINGS)
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
