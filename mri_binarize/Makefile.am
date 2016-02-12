##
## Makefile.am 
##

AM_CFLAGS=-I$(top_srcdir)/include -I$(top_srcdir)/include/dicom
AM_CXXFLAGS=-I$(top_srcdir)/include -I$(top_srcdir)/include/dicom

bin_PROGRAMS = mri_binarize
mri_binarize_SOURCES=mri_binarize.c
mri_binarize_LDADD= $(addprefix $(top_builddir)/, $(LIBS_MGH))
mri_binarize_LDFLAGS=$(OS_LDFLAGS)

# mri_binarize is called by Eugenios hippocampal subfield binaries. 
# Mac OSX systems 10.11 (El Capitan) and greater implemented SIP
# which necessitates a wrapper script to all binaries called by other binaries
# in order to pass the DYLD_LIBRARY_PATH env var.

install-exec-hook:
	cp $(top_builddir)/mri_binarize/mri_binarize $(DESTDIR)$(bindir)/mri_binarize.bin
	echo "#!/bin/bash" > $(DESTDIR)$(bindir)/mri_binarize
	echo "source \$$FREESURFER_HOME/sources.sh" >> $(DESTDIR)$(bindir)/mri_binarize
	echo "mri_binarize.bin \"\$$@\"" >> $(DESTDIR)$(bindir)/mri_binarize

EXCLUDE_FILES=""
include $(top_srcdir)/Makefile.extra
