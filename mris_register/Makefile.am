## 
## Makefile.am 
##
# Makefile.am specify i addtion to the standard

AM_CFLAGS=-I$(top_srcdir)/include -I$(top_srcdir)/include/dicom
AM_LDFLAGS=-L$(top_srcdir)/lib

bin_PROGRAMS = mris_register
mris_register_SOURCES=mris_register.c
mris_register_LDADD= $(addprefix $(top_srcdir)/, $(LIBS_MGH))
mris_register_LDFLAGS=$(OS_LDFLAGS)

