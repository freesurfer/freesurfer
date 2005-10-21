## 
## Makefile.am 
##

AM_CFLAGS=-I$(top_srcdir)/include -I$(top_srcdir)/include/dicom
AM_LDFLAGS=

bin_PROGRAMS = mri_wmfilter
mri_wmfilter_SOURCES=mri_wmfilter.c
mri_wmfilter_LDADD= $(addprefix $(top_srcdir)/, $(LIBS_MGH))
mri_wmfilter_LDFLAGS=$(OS_LDFLAGS)



# Our release target. Include files to be excluded here. They will be
# found and removed after 'make install' is run during the 'make
# release' target.
EXCLUDE_FILES=""
include $(top_srcdir)/Makefile.extra
