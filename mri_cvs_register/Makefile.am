#
# Makefile.am
#

dist_bin_SCRIPTS = mri_cvs_register

foodir=$(DESTDIR)$(bindir)
foo_DATA= \
	elastic_registration_CVS.pl \
	elastic_registration_CVS_hemi.pl \
	getConfig.pm \
	id.xfm \
	mri_cvs_register.settings.txt

EXTRA_DIST=$(foo_DATA)

# Our release target. Include files to be excluded here. They will be
# found and removed after 'make install' is run during the 'make
# release' target.
EXCLUDE_FILES=
include $(top_srcdir)/Makefile.extra
