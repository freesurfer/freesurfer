## 
## Makefile.am 
##
# Makefile.am specify i addtion to the standard

AM_CFLAGS=-I$(top_srcdir)/include -I$(top_srcdir)/include/dicom
AM_LDFLAGS=-L$(top_srcdir)/lib

bin_PROGRAMS = mris_make_template
mris_make_template_SOURCES=mris_make_template.c
mris_make_template_LDADD= $(addprefix $(top_srcdir)/, $(LIBS_MGH))
mris_make_template_LDFLAGS=$(OS_LDFLAGS)

