## 
## Makefile.am 
##

AM_CFLAGS=-I$(top_srcdir)/include -I$(top_srcdir)/include/dicom
AM_LDFLAGS=

bin_PROGRAMS = repair_siemens_file
repair_siemens_file_SOURCES=repair_siemens_file.c
repair_siemens_file_LDADD= $(addprefix $(top_builddir)/, $(LIBS_MGH))
repair_siemens_file_LDFLAGS=$(OS_LDFLAGS)

# Our release target. Include files to be excluded here. They will be
# found and removed after 'make install' is run during the 'make
# release' target.
EXCLUDE_FILES=""
include $(top_srcdir)/Makefile.extra
