##
## Makefile.am 
##

AM_CFLAGS=-I$(top_srcdir)/include
AM_CXXFLAGS=-I$(top_srcdir)/include $(GTS_CFLAGS) -Wno-deprecated

if HAVE_GTS_LIBS
    bin_PROGRAMS = mris_decimate
    mris_decimate_SOURCES=mris_decimate.cpp mris_decimate.h main.cpp
    mris_decimate_LDADD= $(addprefix $(top_builddir)/, $(LIBS_MGH)) $(GTS_LIBS)
    mris_decimate_LDFLAGS=$(OS_LDFLAGS)
endif

# Our release target. Include files to be excluded here. They will be
# found and removed after 'make install' is run during the 'make
# release' target.
EXCLUDE_FILES=
include $(top_srcdir)/Makefile.extra
