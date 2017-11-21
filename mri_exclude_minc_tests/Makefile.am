##
## Makefile.am 
##

AM_CFLAGS=-I$(top_srcdir)/include -I$(top_srcdir)/include/dicom
AM_CXXFLAGS=-I$(top_srcdir)/include -I$(top_srcdir)/include/dicom

bin_PROGRAMS = mri_exclude_minc_tests
mri_exclude_minc_tests_SOURCES=mri_exclude_minc_tests.c
mri_exclude_minc_tests_LDADD= $(addprefix $(top_builddir)/, $(LIBS_MGH))
mri_exclude_minc_tests_LDFLAGS=$(OS_LDFLAGS)

install-exec-hook:
	cp $(top_builddir)/mri_exclude_minc_tests/mri_exclude_minc_tests $(DESTDIR)$(bindir)/mri_exclude_minc_tests.bin
	echo "#!/bin/bash" > $(DESTDIR)$(bindir)/mri_exclude_minc_tests
	echo "source \$$FREESURFER_HOME/sources.sh" >> $(DESTDIR)$(bindir)/mri_exclude_minc_tests
	echo "mri_exclude_minc_tests.bin \"\$$@\"" >> $(DESTDIR)$(bindir)/mri_exclude_minc_tests

EXCLUDE_FILES=""
include $(top_srcdir)/Makefile.extra
