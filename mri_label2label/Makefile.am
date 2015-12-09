## 
## Makefile.am 
##

AM_CFLAGS=-I$(top_srcdir)/include -I$(top_srcdir)/include/dicom
AM_LDFLAGS=

bin_PROGRAMS = mri_label2label
mri_label2label_SOURCES=mri_label2label.c
mri_label2label_LDADD= $(addprefix $(top_builddir)/, $(LIBS_MGH))
mri_label2label_LDFLAGS=$(OS_LDFLAGS)

# mri_label2label is called by qdec. On Mac OSX
# systems 10.11 (El Capitan) and greater implemented SIP
# which necessitates a wrapper script to pass the DYLD_LIBRARY_PATH env var.

install-exec-hook:
	cp $(top_builddir)/mri_label2label/mri_label2label $(DESTDIR)$(bindir)/mri_label2label.bin
	echo "#!/bin/bash" > $(DESTDIR)$(bindir)/mri_label2label
	echo "source \$$FREESURFER_HOME/sources.sh" >> $(DESTDIR)$(bindir)/mri_label2label
	echo "mri_label2label.bin \"\$$@\"" >> $(DESTDIR)$(bindir)/mri_label2label

# Our release target. Include files to be excluded here. They will be
# found and removed after 'make install' is run during the 'make
# release' target.
EXCLUDE_FILES=""
include $(top_srcdir)/Makefile.extra
