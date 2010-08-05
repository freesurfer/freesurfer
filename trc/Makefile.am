##
## Makefile.am 
##

AM_CFLAGS=-I$(top_srcdir)/include
AM_CXXFLAGS=-I$(top_srcdir)/include

bin_PROGRAMS = dmri_spline dmri_track

dmri_spline_SOURCES=dmri_spline.cxx spline.cxx spline.h
dmri_spline_LDADD= $(addprefix $(top_builddir)/, $(LIBS_MGH))
dmri_spline_LDFLAGS=$(OS_LDFLAGS)

dmri_track_SOURCES=dmri_track.cxx coffin.cxx coffin.h \
	bite.cxx bite.h spline.cxx spline.h
dmri_track_LDADD= $(addprefix $(top_builddir)/, $(LIBS_MGH))
dmri_track_LDFLAGS=$(OS_LDFLAGS)

EXCLUDE_FILES=""
include $(top_srcdir)/Makefile.extra
