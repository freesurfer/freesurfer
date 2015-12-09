## 
## Makefile.am 
##

AM_CFLAGS=-I$(top_srcdir)/include -I$(top_srcdir)/include/dicom
AM_LDFLAGS=

bin_PROGRAMS = mri_surfcluster
mri_surfcluster_SOURCES=mri_surfcluster.c
mri_surfcluster_LDADD= $(addprefix $(top_builddir)/, $(LIBS_MGH))
mri_surfcluster_LDFLAGS=$(OS_LDFLAGS)

# mri_surfcluster is called by qdec. On Mac OSX
# systems 10.11 (El Capitan) and greater implemented SIP
# which necessitates a wrapper script to pass the DYLD_LIBRARY_PATH env var.

install-exec-hook:
	cp $(top_builddir)/mri_surfcluster/mri_surfcluster $(DESTDIR)$(bindir)/mri_surfcluster.bin
	echo "#!/bin/bash" > $(DESTDIR)$(bindir)/mri_surfcluster
	echo "source \$$FREESURFER_HOME/sources.sh" >> $(DESTDIR)$(bindir)/mri_surfcluster
	echo "mri_surfcluster.bin \"\$$@\"" >> $(DESTDIR)$(bindir)/mri_surfcluster

# Our release target. Include files to be excluded here. They will be
# found and removed after 'make install' is run during the 'make
# release' target.
EXCLUDE_FILES=""
include $(top_srcdir)/Makefile.extra
