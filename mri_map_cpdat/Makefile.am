##
## Makefile.am 
##

AM_CFLAGS=-I$(top_srcdir)/include
AM_CXXFLAGS=-I$(top_srcdir)/include

bin_PROGRAMS = mri_map_cpdat
mri_map_cpdat_SOURCES=mri_map_cpdat.cpp
mri_map_cpdat_LDADD= $(addprefix $(top_builddir)/, $(LIBS_MGH))
mri_map_cpdat_LDFLAGS=$(OS_LDFLAGS)


# Our release target. Include files to be excluded here. They will be
# found and removed after 'make install' is run during the 'make
# release' target.
EXCLUDE_FILES=
include $(top_srcdir)/Makefile.extra
