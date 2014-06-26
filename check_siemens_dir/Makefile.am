## 
# Makefile.am
##

AM_CFLAGS=-I$(top_srcdir)/include -I$(top_srcdir)/include/dicom
AM_LDFLAGS=

bin_PROGRAMS = check_siemens_dir
check_siemens_dir_SOURCES=check_siemens_dir.c
check_siemens_dir_LDADD= $(addprefix $(top_builddir)/, $(LIBS_MGH))
check_siemens_dir_LDFLAGS=$(OS_LDFLAGS)

# Our release target. Include files to be excluded here. They will be
# found and removed after 'make install' is run during the 'make
# release' target.
EXCLUDE_FILES=""
include $(top_srcdir)/Makefile.extra
