## 
## Makefile.am 
##
# Makefile.am specify i addtion to the standard

AM_CFLAGS=-I$(top_srcdir)/include -I$(top_srcdir)/include/dicom
AM_LDFLAGS=-L$(top_srcdir)/lib

bin_PROGRAMS = mris_w_to_curv
mris_w_to_curv_SOURCES=mris_w_to_curv.c
mris_w_to_curv_LDADD= $(addprefix $(top_srcdir)/, $(LIBS_MGH))
mris_w_to_curv_LDFLAGS=$(OS_LDFLAGS)

