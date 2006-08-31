##
## Makefile.am 
##

AM_CFLAGS=-I$(top_srcdir)/include -I$(top_srcdir)/include/dicom $(GSL_CFLAGS)
AM_CXXFLAGS=-I$(top_srcdir)/include -I$(top_srcdir)/include/dicom

bin_PROGRAMS = mris_niters2fwhm
mris_niters2fwhm_SOURCES=mris_niters2fwhm.c
mris_niters2fwhm_LDADD= $(addprefix $(top_srcdir)/, $(LIBS_MGH))
mris_niters2fwhm_LDFLAGS=$(OS_LDFLAGS)

# Our release target. Include files to be excluded here. They will be
# found and removed after 'make install' is run during the 'make
# release' target.
EXCLUDE_FILES=""
include $(top_srcdir)/Makefile.extra
