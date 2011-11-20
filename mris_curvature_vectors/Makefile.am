## 
## Makefile.am 
##
# Makefile.am specify i addtion to the standard

AM_CFLAGS=-I$(top_srcdir)/include -I$(top_srcdir)/include/dicom
AM_LDFLAGS=

bin_PROGRAMS = mris_curvature_vectors
mris_curvature_vectors_SOURCES=mris_curvature_vectors.c
mris_curvature_vectors_LDADD= $(addprefix $(top_builddir)/, $(LIBS_MGH))
mris_curvature_vectors_LDFLAGS=$(OS_LDFLAGS)



# Our release target. Include files to be excluded here. They will be
# found and removed after 'make install' is run during the 'make
# release' target.
EXCLUDE_FILES=""
include $(top_srcdir)/Makefile.extra
