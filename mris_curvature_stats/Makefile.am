## 
## Makefile.am 
##

AM_CFLAGS=-I$(top_srcdir)/include -I$(top_srcdir)/include/dicom
AM_LDFLAGS=

bin_PROGRAMS = mris_curvature_stats
mris_curvature_stats_SOURCES=mris_curvature_stats.c
mris_curvature_stats_LDADD= $(addprefix $(top_builddir)/, $(LIBS_MGH))
mris_curvature_stats_LDFLAGS=$(OS_LDFLAGS)

# Our release target. Include files to be excluded here. They will be
# found and removed after 'make install' is run during the 'make
# release' target.
EXCLUDE_FILES=""
include $(top_srcdir)/Makefile.extra
