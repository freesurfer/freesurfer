#
# Makefile.am
#

dist_bin_SCRIPTS = mri_cvs_register

foodir=$(DESTDIR)$(bindir)
foo_DATA= \
	elastic_registration_CVS.pl \
	getConfig.pm \
	id.xfm \
	set-mpi17-3.txt

# Our release target. Include files to be excluded here. They will be
# found and removed after 'make install' is run during the 'make
# release' target.
EXCLUDE_FILES=
include $(top_srcdir)/Makefile.extra
