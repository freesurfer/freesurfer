## 
## Makefile.am 
##
AM_CFLAGS=-I$(top_srcdir)/include -I$(top_srcdir)/include/dicom
AM_CXXFLAGS=-I$(top_srcdir)/include -I$(top_srcdir)/include/dicom

bin_PROGRAMS = dngtester
EXTRASRC=
dngtester_SOURCES=dngtester.c  $(EXTRASRC)
dngtester_LDADD= $(addprefix $(top_srcdir)/, $(LIBS_MGH))
dngtester_LDFLAGS=$(OS_LDFLAGS)
