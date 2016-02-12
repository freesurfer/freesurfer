##
## Makefile.am 
##

AM_CFLAGS=-I$(top_srcdir)/include

bin_PROGRAMS = mri_concatenate_lta
mri_concatenate_lta_SOURCES=mri_concatenate_lta.c
mri_concatenate_lta_LDADD= $(addprefix $(top_builddir)/, $(LIBS_MGH))
mri_concatenate_lta_LDFLAGS=$(OS_LDFLAGS)

# trick to get test data into the build directory during make distcheck
foo_SOURCES=
foo:
	if ! test -f "$(top_builddir)/mri_concatenate_lta/testdata.tar.gz"; then \
		cp -vf $(top_srcdir)/mri_concatenate_lta/testdata.tar.gz \
					$(top_builddir)/mri_concatenate_lta/; fi;

check_PROGRAMS=foo

TESTS=test_mri_concatenate_lta

EXTRA_DIST=test_mri_concatenate_lta testdata.tar.gz

# mri_concatenate_lta is called by Eugenios hippocampal subfield binaries. 
# Mac OSX systems 10.11 (El Capitan) and greater implemented SIP
# which necessitates a wrapper script to all binaries called by other binaries
# in order to pass the DYLD_LIBRARY_PATH env var.

install-exec-hook:
	cp $(top_builddir)/mri_concatenate_lta/mri_concatenate_lta $(DESTDIR)$(bindir)/mri_concatenate_lta.bin
	echo "#!/bin/bash" > $(DESTDIR)$(bindir)/mri_concatenate_lta
	echo "source \$$FREESURFER_HOME/sources.sh" >> $(DESTDIR)$(bindir)/mri_concatenate_lta
	echo "mri_concatenate_lta.bin \"\$$@\"" >> $(DESTDIR)$(bindir)/mri_concatenate_lta

# Our release target. Include files to be excluded here. They will be
# found and removed after 'make install' is run during the 'make
# release' target.
EXCLUDE_FILES=""
include $(top_srcdir)/Makefile.extra
