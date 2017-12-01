##
## Makefile.am 
##

AM_CFLAGS=-I$(top_srcdir)/include -I$(top_srcdir)/include/dicom
AM_CXXFLAGS=-I$(top_srcdir)/include -I$(top_srcdir)/include/dicom

bin_PROGRAMS = test_exclude_minc
test_exclude_minc_SOURCES=test_exclude_minc.c
test_exclude_minc_LDADD= $(addprefix $(top_builddir)/, $(LIBS_MGH))
test_exclude_minc_LDFLAGS=$(OS_LDFLAGS)

install-exec-hook:
	cp $(top_builddir)/test_exclude_minc/test_exclude_minc $(DESTDIR)$(bindir)/test_exclude_minc.bin
	echo "#!/bin/bash" > $(DESTDIR)$(bindir)/test_exclude_minc
	echo "source \$$FREESURFER_HOME/sources.sh" >> $(DESTDIR)$(bindir)/test_exclude_minc
	echo "test_exclude_minc.bin \"\$$@\"" >> $(DESTDIR)$(bindir)/test_exclude_minc

EXCLUDE_FILES=""
include $(top_srcdir)/Makefile.extra
