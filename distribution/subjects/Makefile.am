#
# Makefile.am
#

foodir=$(prefix)/subjects
foo_datafiles=README sample-001.mgz sample-002.mgz

install-data-local:
	test -z $(DESTDIR)$(foodir) || $(mkdir_p) $(DESTDIR)$(foodir)
		$(INSTALL_DATA) $(foo_datafiles) $(DESTDIR)$(foodir)

# Our release target. Include files to be excluded here. They will be
# found and removed after 'make install' is run during the 'make
# release' target.
EXCLUDE_FILES=
include $(top_srcdir)/Makefile.extra
