## 
## Makefile.am 
##
# Makefile.am specify i addtion to the standard

AM_CFLAGS=-I$(top_srcdir)/include -I$(top_srcdir)/include/dicom
AM_CXXFLAGS=-I$(top_srcdir)/include -I$(top_srcdir)/include/dicom
AM_LDFLAGS=-L$(top_srcdir)/lib

bin_PROGRAMS = mris_distance_to_label
mris_distance_to_label_SOURCES=mris_distance_to_label.cpp mri_fastmarching.cpp
mris_distance_to_label_LDADD= $(addprefix $(top_srcdir)/, $(LIBS_MGH))
mris_distance_to_label_LDFLAGS=$(OS_LDFLAGS)

