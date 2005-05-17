## 
## Makefile.am 
AM_CFLAGS=-I$(top_srcdir)/include -I$(top_srcdir)/include/dicom
AM_CXXFLAGS=-I$(top_srcdir)/include -I$(top_srcdir)/include/dicom

bin_PROGRAMS = label_area
label_area_SOURCES=label_area.c
label_area_LDADD= $(addprefix $(top_srcdir)/, $(LIBS_MGH))
label_area_LDFLAGS=$(OS_LDFLAGS)

check_PROGRAMS= myown

myown_SOURCES=myown.c
myown_LDADD=$(addprefix $(top_srcdir)/, $(LIBS_MGH))



