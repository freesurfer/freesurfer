## 
## Makefile.am 
##
# Makefile.am specify i addtion to the standard
AM_CFLAGS=-I$(top_srcdir)/include -I$(top_srcdir)/include/dicom
AM_LDFLAGS=

bin_PROGRAMS = mris_fix_topology mris_remove_intersection
mris_fix_topology_SOURCES=mris_fix_topology.c
mris_fix_topology_LDADD= $(addprefix $(top_srcdir)/, $(LIBS_MGH))
mris_fix_topology_LDFLAGS=$(OS_LDFLAGS)

mris_remove_intersection_SOURCES=mris_remove_intersection.c
mris_remove_intersection_LDADD= $(addprefix $(top_srcdir)/, $(LIBS_MGH))
mris_remove_intersection_LDFLAGS=$(OS_LDFLAGS)

# Our release target. Include files to be excluded here. They will be
# found and removed after 'make install' is run during the 'make
# release' target.
EXCLUDE_FILES=""
include $(top_srcdir)/Makefile.extra
