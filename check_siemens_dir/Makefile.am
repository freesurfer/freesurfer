## 
# Makefile.am - specify in addition to the standard
##
AM_CFLAGS=-I$(top_srcdir)/include \
	-I$(top_srcdir)/include/dicom
AM_LDFLAGS=-L$(top_srcdir)/lib

bin_PROGRAMS = check_siemens_dir
check_siemens_dir_SOURCES=check_siemens_dir.c
check_siemens_dir_LDADD= $(addprefix $(top_srcdir)/, $(LIBS_MGH))
check_siemens_dir_LDFLAGS=$(OS_LDFLAGS)

