## 
## Makefile.am 
##

AM_CFLAGS=-I$(top_srcdir)/include -I$(top_srcdir)/include/dicom
AM_CXXFLAGS=-I$(top_srcdir)/include -I$(top_srcdir)/include/dicom

bin_PROGRAMS = mri_concat
mri_concat_SOURCES=mri_concat.c
mri_concat_LDADD= $(addprefix $(top_builddir)/, $(LIBS_MGH))
mri_concat_LDFLAGS=$(OS_LDFLAGS)

# mri_concat is called by qdec. On Mac OSX
# systems 10.11 (El Capitan) and greater implemented SIP
# which necessitates a wrapper script to pass the DYLD_LIBRARY_PATH env var.

install-exec-hook:
	cp $(top_builddir)/mri_concat/mri_concat $(DESTDIR)$(bindir)/mri_concat.bin
	echo "#!/bin/bash" > $(DESTDIR)$(bindir)/mri_concat
	echo "source \$$FREESURFER_HOME/sources.sh" >> $(DESTDIR)$(bindir)/mri_concat
	echo "mri_concat.bin \"\$$@\"" >> $(DESTDIR)$(bindir)/mri_concat

# Our release target. Include files to be excluded here. They will be
# found and removed after 'make install' is run during the 'make
# release' target.
EXCLUDE_FILES=""
include $(top_srcdir)/Makefile.extra
