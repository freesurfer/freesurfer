## 
## Makefile.am 
##
# Makefile.am specify i addtion to the standard

AM_CFLAGS=-I$(top_srcdir)/include -I$(top_srcdir)/include/dicom
AM_LDFLAGS=-L$(top_srcdir)/lib

bin_PROGRAMS = mris_twoclass
mris_twoclass_SOURCES=mris_twoclass.c
mris_twoclass_LDADD= $(addprefix $(top_srcdir)/, $(LIBS_MGH))
mris_twoclass_LDFLAGS=$(OS_LDFLAGS)

