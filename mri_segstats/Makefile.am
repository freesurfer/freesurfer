## 
## Makefile.am 
##

AM_CFLAGS=-I$(top_srcdir)/include -I$(top_srcdir)/include/dicom
AM_CXXFLAGS=-I$(top_srcdir)/include -I$(top_srcdir)/include/dicom

bin_PROGRAMS = mri_segstats
mri_segstats_SOURCES=mri_segstats.c
mri_segstats_LDADD= $(addprefix $(top_srcdir)/, $(LIBS_MGH))
mri_segstats_LDFLAGS=$(OS_LDFLAGS)

check_PROGRAMS= myown

myown_SOURCES=myown.c
myown_LDADD=$(addprefix $(top_srcdir)/, $(LIBS_MGH))
