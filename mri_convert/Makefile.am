## 
## Makefile.am 
##
# Makefile.am specify i addtion to the standard
AM_CPPFLAGS=-I$(top_srcdir)/include -I$(top_srcdir)/include/dicom
AM_LDFLAGS=

bin_PROGRAMS = mri_convert
mri_convert_SOURCES=mri_convert.c
mri_convert_LDADD=
mri_convert_DEPENDENCIES=../utils/libutils.a
mri_convert_LDFLAGS=-L$(top_srcdir)/lib/$(OS) -lutils -lrgb -ldicom -lunix -lhipsstubs $(LIBS)
