## 
## Makefile.am 
##
# Makefile.am specify i addtion to the standard

AM_CFLAGS=-I$(top_srcdir)/include -I$(top_srcdir)/include/dicom
AM_LDFLAGS=-L$(top_srcdir)/lib

bin_PROGRAMS = mris_jacobian
mris_jacobian_SOURCES=mris_jacobian.c
mris_jacobian_LDADD= $(addprefix $(top_srcdir)/, $(LIBS_MGH))
mris_jacobian_LDFLAGS=$(OS_LDFLAGS)

