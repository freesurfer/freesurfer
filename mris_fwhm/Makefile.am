##
## Makefile.am 
##

AM_CFLAGS=-I$(top_srcdir)/include
AM_CXXFLAGS=-I$(top_srcdir)/include

bin_PROGRAMS = mris_fwhm
mris_fwhm_SOURCES=mris_fwhm.c
mris_fwhm_LDADD= $(addprefix $(top_srcdir)/, $(LIBS_MGH))
mris_fwhm_LDFLAGS=$(OS_LDFLAGS)

# Our release target. Include files to be excluded here. They will be
# found and removed after 'make install' is run during the 'make
# release' target.
EXCLUDE_FILES=
include $(top_srcdir)/Makefile.extra
