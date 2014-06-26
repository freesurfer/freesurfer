##
## Makefile.am 
##

AM_CFLAGS=-I$(top_srcdir)/include -I$(top_srcdir)/include/dicom
AM_CXXFLAGS=-I$(top_srcdir)/include -I$(top_srcdir)/include/dicom

bin_PROGRAMS = mri_z2p
mri_z2p_SOURCES=mri_z2p.c
mri_z2p_LDADD= $(addprefix $(top_builddir)/, $(LIBS_MGH))
mri_z2p_LDFLAGS=$(OS_LDFLAGS)

check_PROGRAMS=

# Our release target. Include files to be excluded here. They will be
# found and removed after 'make install' is run during the 'make
# release' target.
EXCLUDE_FILES=""
include $(top_srcdir)/Makefile.extra
