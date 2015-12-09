## 
## Makefile.am 
##

AM_CFLAGS=-I$(top_srcdir)/include
AM_CXXFLAGS=-I$(top_srcdir)/include

bin_PROGRAMS = mri_glmfit
mri_glmfit_SOURCES=mri_glmfit.c $(top_srcdir)/fsgdf/fsgdf.c
mri_glmfit_LDADD=$(addprefix $(top_builddir)/, $(LIBS_MGH))
mri_glmfit_LDFLAGS=$(OS_LDFLAGS)

# trick to get test data into the build directory
foo_SOURCES=
foo:
	if ! test -f "$(top_builddir)/mri_glmfit/teststats.tar.gz"; then \
		cp -v $(top_srcdir)/mri_glmfit/teststats.tar.gz \
					$(top_builddir)/mri_glmfit/; fi;

check_PROGRAMS=foo

TESTS=test_mri_glmfit

# mri_glmfit is called by qdec. On Mac OSX
# systems 10.11 (El Capitan) and greater implemented SIP
# which necessitates a wrapper script to pass the DYLD_LIBRARY_PATH env var.

install-exec-hook:
	cp $(top_builddir)/mri_glmfit/mri_glmfit $(DESTDIR)$(bindir)/mri_glmfit.bin
	echo "#!/bin/bash" > $(DESTDIR)$(bindir)/mri_glmfit
	echo "source \$$FREESURFER_HOME/sources.sh" >> $(DESTDIR)$(bindir)/mri_glmfit
	echo "mri_glmfit.bin \"\$$@\"" >> $(DESTDIR)$(bindir)/mri_glmfit

EXTRA_DIST=test_mri_glmfit teststats.tar.gz

# Our release target. Include files to be excluded here. They will be
# found and removed after 'make install' is run during the 'make
# release' target.
EXCLUDE_FILES=
include $(top_srcdir)/Makefile.extra
