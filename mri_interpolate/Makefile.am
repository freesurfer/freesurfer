##
## Makefile.am 
##

AM_CFLAGS=-I$(top_srcdir)/include
AM_CXXFLAGS=-I$(top_srcdir)/include

bin_PROGRAMS = mri_interpolate
mri_interpolate_SOURCES=mri_interpolate.c
mri_interpolate_LDADD= $(addprefix $(top_srcdir)/, $(LIBS_MGH))
mri_interpolate_LDFLAGS=$(OS_LDFLAGS)

# Our release target. Include files to be excluded here. They will be
# found and removed after 'make install' is run during the 'make
# release' target.
EXCLUDE_FILES=
include $(top_srcdir)/Makefile.extra
