## 
## Makefile.am 
##
# Makefile.am specify i addtion to the standard

AM_CFLAGS=-I$(top_srcdir)/include -I$(top_srcdir)/include/dicom
AM_LDFLAGS=-L$(top_srcdir)/lib

bin_PROGRAMS = repair_siemens_file
repair_siemens_file_SOURCES=repair_siemens_file.c
repair_siemens_file_LDADD= $(addprefix $(top_srcdir)/, $(LIBS_MGH))
repair_siemens_file_LDFLAGS=$(OS_LDFLAGS)

