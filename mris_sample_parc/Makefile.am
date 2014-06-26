## 
## Makefile.am 
##

AM_CFLAGS=-I$(top_srcdir)/include -I$(top_srcdir)/include/dicom
AM_LDFLAGS=

bin_PROGRAMS = mris_sample_parc
mris_sample_parc_SOURCES=mris_sample_parc.c
mris_sample_parc_LDADD= $(addprefix $(top_builddir)/, $(LIBS_MGH))
mris_sample_parc_LDFLAGS=$(OS_LDFLAGS)

# Our release target. Include files to be excluded here. They will be
# found and removed after 'make install' is run during the 'make
# release' target.
EXCLUDE_FILES=""
include $(top_srcdir)/Makefile.extra
