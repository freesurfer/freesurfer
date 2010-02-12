## 
## Makefile.am 
##

AM_CFLAGS=-I$(top_srcdir)/include
AM_CXXFLAGS=-I$(top_srcdir)/include

bin_PROGRAMS = mri_mi_nl_align
mri_mi_nl_align_SOURCES=mri_mi_nl_align.c gcam_mi.cpp
mri_mi_nl_align_LDADD= $(addprefix $(top_builddir)/, $(LIBS_MGH))
mri_mi_nl_align_LDFLAGS=$(OS_LDFLAGS)

EXTRA_DIST=gcam_mi.h

# Our release target. Include files to be excluded here. They will be
# found and removed after 'make install' is run during the 'make
# release' target.
EXCLUDE_FILES=""
include $(top_srcdir)/Makefile.extra
