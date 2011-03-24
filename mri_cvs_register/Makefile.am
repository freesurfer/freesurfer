#
# Makefile.am
#

dist_bin_SCRIPTS = mri_cvs_register mri_cvs_check mri_cvs_data_copy

foodir=$(prefix)/bin
foo_DATA= \
	elastic_registration_CVS.pl \
	elastic_registration_CVS_hemi.pl \
	getConfig.pm \
	id.xfm \
	mri_cvs_register.settings.txt \
	mri_cvs_requiredfiles.txt

xmldir=$(prefix)/docs/xml
xml_DATA=mri_cvs_register.help.xml

install-exec-hook:
	$(top_builddir)/utils/fsPrintHelp $(srcdir)/mri_cvs_register.help.xml \
		>> $(DESTDIR)$(bindir)/mri_cvs_register

TESTS =  $(top_builddir)/scripts/help_xml_validate

EXTRA_DIST=$(foo_DATA)

# Our release target. Include files to be excluded here. They will be
# found and removed after 'make install' is run during the 'make
# release' target.
EXCLUDE_FILES=
include $(top_srcdir)/Makefile.extra
