## 
## Makefile.am 
##
# Makefile.am specify i addtion to the standard

AM_CXXFLAGS=-I$(top_srcdir)/include -I$(top_srcdir)/include/dicom
AM_LDFLAGS=-L$(top_srcdir)/lib

bin_PROGRAMS = regdat2xfm
regdat2xfm_SOURCES=regdat2xfm.cpp
regdat2xfm_LDADD= $(addprefix $(top_srcdir)/, $(LIBS_MGH))
regdat2xfm_LDFLAGS=$(OS_LDFLAGS)

