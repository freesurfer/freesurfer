## 
## Makefile.am 
##
# Makefile.am specify i addtion to the standard

AM_CFLAGS=-I$(top_srcdir)/include -I$(top_srcdir)/include/dicom
AM_LDFLAGS=-L$(top_srcdir)/lib

bin_PROGRAMS = rbftest
rbftest_SOURCES=rbftest.c
rbftest_LDADD= $(addprefix $(top_srcdir)/, $(LIBS_MGH))
rbftest_LDFLAGS=$(OS_LDFLAGS)

