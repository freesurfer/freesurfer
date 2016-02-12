## 
## Makefile.am 
##

AM_CFLAGS=-I$(top_srcdir)/include
AM_LDFLAGS=

bin_PROGRAMS = mri_mask

mri_mask_SOURCES=mri_mask.c mri_mask.help.xml.h
mri_mask_LDADD= $(addprefix $(top_builddir)/, $(LIBS_MGH))
mri_mask_LDFLAGS=$(OS_LDFLAGS)

foodir=$(prefix)/docs/xml
foo_DATA=mri_mask.help.xml
foo2dir=$(prefix)/docs/html
foo2_DATA=mri_mask.help.xml.html

TESTS=$(top_builddir)/scripts/help_xml_validate

clean-local:
	rm -f $(BUILT_SOURCES)

EXTRA_DIST=$(foo_DATA) $(BUILT_SOURCES)

# mri_mask is called by Eugenios hippocampal subfield binaries. 
# Mac OSX systems 10.11 (El Capitan) and greater implemented SIP
# which necessitates a wrapper script to all binaries called by other binaries
# in order to pass the DYLD_LIBRARY_PATH env var.

install-exec-hook:
	cp $(top_builddir)/mri_mask/mri_mask $(DESTDIR)$(bindir)/mri_mask.bin
	echo "#!/bin/bash" > $(DESTDIR)$(bindir)/mri_mask
	echo "source \$$FREESURFER_HOME/sources.sh" >> $(DESTDIR)$(bindir)/mri_mask
	echo "mri_mask.bin \"\$$@\"" >> $(DESTDIR)$(bindir)/mri_mask

# Our release target. Include files to be excluded here. They will be
# found and removed after 'make install' is run during the 'make
# release' target.
EXCLUDE_FILES=""
include $(top_srcdir)/Makefile.extra

BUILT_SOURCES=mri_mask.help.xml.h mri_mask.help.xml.html

