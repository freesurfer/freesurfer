## 
## Makefile.am 
##
AM_CFLAGS=-I$(top_srcdir)/include -I$(top_srcdir)/include/dicom
AM_CXXFLAGS=-I$(top_srcdir)/include -I$(top_srcdir)/include/dicom

bin_PROGRAMS = mris_volume
mris_volume_SOURCES=mris_volume.c
mris_volume_LDADD= $(addprefix $(top_srcdir)/, $(LIBS_MGH))
mris_volume_LDFLAGS=$(OS_LDFLAGS)


