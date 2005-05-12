## 
## Makefile.am 
##
# Makefile.am specify i addtion to the standard

AM_CFLAGS=-I$(top_srcdir)/include -I$(top_srcdir)/include/dicom
AM_CXXFLAGS=-I$(top_srcdir)/include -I$(top_srcdir)/include/dicom

bin_PROGRAMS = mri_distance_transform
mri_distance_transform_SOURCES=mri_distance_transform.cpp
mri_distance_transform_LDADD= $(addprefix $(top_srcdir)/, $(LIBS_MGH) utilscpp/libutilscpp.a )
mri_distance_transform_LDFLAGS=$(OS_LDFLAGS)

