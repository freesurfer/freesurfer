## 
## Makefile.am 
##
# Makefile.am specify i addtion to the standard

AM_CFLAGS=-I$(top_srcdir)/include -I$(top_srcdir)/include/dicom
AM_LDFLAGS=-L$(top_srcdir)/lib

bin_PROGRAMS=histo_compute_joint_density
histo_compute_joint_density_SOURCES=histo_compute_joint_density.c
histo_compute_joint_density_LDADD= $(addprefix $(top_srcdir)/, $(LIBS_MGH))
histo_compute_joint_density_LDFLAGS=$(OS_LDFLAGS)


