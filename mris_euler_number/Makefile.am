## 
## Makefile.am 
##
# Makefile.am specify i addtion to the standard

AM_CFLAGS=-I$(top_srcdir)/include -I$(top_srcdir)/include/dicom
AM_LDFLAGS=-L$(top_srcdir)/lib

bin_PROGRAMS = mris_euler_number
mris_euler_number_SOURCES=mris_euler_number.c
mris_euler_number_LDADD= $(addprefix $(top_srcdir)/, $(LIBS_MGH))
mris_euler_number_LDFLAGS=$(OS_LDFLAGS)

