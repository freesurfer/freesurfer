#
# Makefile.am
#

dist_bin_SCRIPTS = mri_cvs_register mri_cvs_check mri_cvs_data_copy

foodir=$(bindir)
foo_DATA= \
	id.xfm \
	mri_cvs_requiredfiles.txt

xmldir=$(prefix)/docs/xml
xml_DATA=mri_cvs_register.help.xml
htmldir=$(prefix)/docs/html
html_DATA=mri_cvs_register.help.xml.html

install-exec-hook:
	$(top_builddir)/utils/fsPrintHelp $(srcdir)/mri_cvs_register.help.xml \
		>> $(DESTDIR)$(bindir)/mri_cvs_register

TESTS =  $(top_builddir)/scripts/help_xml_validate

EXTRA_DIST=$(foo_DATA) $(xml_DATA)

# Our release target. Include files to be excluded here. They will be
# found and removed after 'make install' is run during the 'make
# release' target.
EXCLUDE_FILES=
include $(top_srcdir)/Makefile.extra

clean-local:
	rm -f *.help.xml.html
