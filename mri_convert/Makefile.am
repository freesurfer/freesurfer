## 
## Makefile.am 
##

AM_CFLAGS=-I$(top_srcdir)/include -I$(top_srcdir)/include/dicom

bin_PROGRAMS = mri_convert mri_make_uchar
mri_convert_SOURCES=mri_convert.c mri_convert.help.xml.h
mri_convert_LDADD= $(addprefix $(top_builddir)/, $(LIBS_MGH))
mri_convert_LDFLAGS=$(OS_LDFLAGS)

mri_make_uchar_SOURCES=mri_make_uchar.c
mri_make_uchar_LDADD= $(addprefix $(top_builddir)/, $(LIBS_MGH))
mri_make_uchar_LDFLAGS=$(OS_LDFLAGS)

foodir=$(prefix)/docs/xml
foo_DATA=mri_convert.help.xml
foo2dir=$(prefix)/docs/html
foo2_DATA=mri_convert.help.xml.html

# trick to get test data into the build directory
foo_SOURCES=
foo:
	if ! test -f "$(top_builddir)/mri_convert/testdata.tar.gz"; then \
		cp -v $(top_srcdir)/mri_convert/testdata.tar.gz \
					$(top_builddir)/mri_convert/; fi;

check_PROGRAMS=foo

TESTS=test_mri_convert $(top_builddir)/scripts/help_xml_validate

clean-local:
	rm -f $(BUILT_SOURCES)

EXTRA_DIST=test_mri_convert testdata.tar.gz $(foo_DATA) $(BUILT_SOURCES)

# mri_convert is called by Eugenio's hippocampal subfield code.
# Mac OSX systems 10.11 (El Capitan) and greater implemented SIP
# which necessitates a wrapper script for all binaries called by 
# other binaries to pass the DYLD_LIBRARY_PATH env var.

install-exec-hook:
	cp $(top_builddir)/mri_convert/mri_convert $(DESTDIR)$(bindir)/mri_convert.bin
	echo "#!/bin/bash" > $(DESTDIR)$(bindir)/mri_convert
	echo "source \$$FREESURFER_HOME/sources.sh" >> $(DESTDIR)$(bindir)/mri_convert
	echo "mri_convert.bin \"\$$@\"" >> $(DESTDIR)$(bindir)/mri_convert

# Our release target. Include files to be excluded here. They will be
# found and removed after 'make install' is run during the 'make
# release' target.
EXCLUDE_FILES=
include $(top_srcdir)/Makefile.extra

BUILT_SOURCES=mri_convert.help.xml.h mri_convert.help.xml.html
