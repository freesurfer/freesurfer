## 
## Makefile.am 
##
# Makefile.am specify i addtion to the standard

AM_CFLAGS=-I$(top_srcdir)/include -I$(top_srcdir)/include/dicom
AM_LDFLAGS=

bin_PROGRAMS = mri_ms_fitparms
mri_ms_fitparms_SOURCES=mri_ms_fitparms.c
mri_ms_fitparms_LDADD= $(addprefix $(top_srcdir)/, $(LIBS_MGH))
mri_ms_fitparms_LDFLAGS=$(OS_LDFLAGS)



# Our release target. Include files to be excluded here. They will be
# found and removed after 'make install' is run during the 'make
# release' target.
EXCLUDE_FILES=""
include $(top_srcdir)/Makefile.extra
