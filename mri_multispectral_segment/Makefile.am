##
## Makefile.am 
##

AM_CFLAGS=-I$(top_srcdir)/include -I$(top_srcdir)/include/dicom
AM_CXXFLAGS=-I$(top_srcdir)/include -I$(top_srcdir)/include/dicom

bin_PROGRAMS = mri_multispectral_segment
mri_multispectral_segment_SOURCES=mri_multispectral_segment.c
mri_multispectral_segment_LDADD= $(addprefix $(top_srcdir)/, $(LIBS_MGH))
mri_multispectral_segment_LDFLAGS=$(OS_LDFLAGS)

# Our release target. Include files to be excluded here. They will be
# found and removed after 'make install' is run during the 'make
# release' target.
EXCLUDE_FILES=""
include $(top_srcdir)/Makefile.extra
