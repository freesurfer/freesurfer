## 
## Makefile.am 
##

AM_CFLAGS=-I$(top_srcdir)/include -I$(top_srcdir)/include/dicom
AM_LDFLAGS=

bin_PROGRAMS = mris_compute_acorr
mris_compute_acorr_SOURCES=mris_compute_acorr.c
mris_compute_acorr_LDADD= $(addprefix $(top_builddir)/, $(LIBS_MGH))
mris_compute_acorr_LDFLAGS=$(OS_LDFLAGS)

# Our release target. Include files to be excluded here. They will be
# found and removed after 'make install' is run during the 'make
# release' target.
EXCLUDE_FILES=""
include $(top_srcdir)/Makefile.extra
