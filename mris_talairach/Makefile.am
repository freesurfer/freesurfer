## 
## Makefile.am 
##
# Makefile.am specify i addtion to the standard

AM_CFLAGS=-I$(top_srcdir)/include -I$(top_srcdir)/include/dicom
AM_LDFLAGS=-L$(top_srcdir)/lib

bin_PROGRAMS = mris_talairach
mris_talairach_SOURCES=mris_talairach.c
mris_talairach_LDADD= $(addprefix $(top_srcdir)/, $(LIBS_MGH))
mris_talairach_LDFLAGS=$(OS_LDFLAGS)

