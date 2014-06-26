##
## Makefile.am 
##

AM_CFLAGS=-I$(top_srcdir)/include
AM_CXXFLAGS=-I$(top_srcdir)/include

bin_PROGRAMS = mri_multiscale_segment
mri_multiscale_segment_SOURCES=mri_multiscale_segment.c
mri_multiscale_segment_LDADD= $(addprefix $(top_builddir)/, $(LIBS_MGH))
mri_multiscale_segment_LDFLAGS=$(OS_LDFLAGS)

# Our release target. Include files to be excluded here. They will be
# found and removed after 'make install' is run during the 'make
# release' target.
EXCLUDE_FILES=
include $(top_srcdir)/Makefile.extra
