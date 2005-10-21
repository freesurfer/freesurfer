## 
## Makefile.am 
##
# Makefile.am specify i addtion to the standard

AM_CFLAGS=-I$(top_srcdir)/include -I$(top_srcdir)/include/dicom
AM_LDFLAGS=

bin_PROGRAMS = mri_ca_normalize
mri_ca_normalize_SOURCES=mri_ca_normalize.c
mri_ca_normalize_LDADD= $(addprefix $(top_srcdir)/, $(LIBS_MGH))
mri_ca_normalize_LDFLAGS=$(OS_LDFLAGS)


# Our release target. Include files to be excluded here. They will be
# found and removed after 'make install' is run during the 'make
# release' target.
EXCLUDE_FILES=""
include $(top_srcdir)/Makefile.extra
