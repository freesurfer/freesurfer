## 
## Makefile.am 
##
# Makefile.am specify i addtion to the standard

AM_CFLAGS=-I$(top_srcdir)/include -I$(top_srcdir)/include/dicom
AM_LDFLAGS=-L$(top_srcdir)/lib

bin_PROGRAMS = mris_fix_topology
mris_fix_topology_SOURCES=mris_fix_topology.c
mris_fix_topology_LDADD= $(addprefix $(top_srcdir)/, $(LIBS_MGH))
#mris_fix_topology_LDFLAGS=$(OS_LDFLAGS)

