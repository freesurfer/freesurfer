## 
## Makefile.am - specify i addition to the standard
##

AM_CFLAGS=-I$(top_srcdir)/include -I$(top_srcdir)/include/dicom
AM_LDFLAGS=

bin_PROGRAMS = mkxsubjreg
mkxsubjreg_SOURCES=mkxsubjreg.c
mkxsubjreg_LDADD= $(addprefix $(top_srcdir)/, $(LIBS_MGH))
mkxsubjreg_LDFLAGS=$(OS_LDFLAGS)



# Our release target. Include files to be excluded here. They will be
# found and removed after 'make install' is run during the 'make
# release' target.
EXCLUDE_FILES=""
include $(top_srcdir)/Makefile.extra
