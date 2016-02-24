#
# Makefile.am
#

foodir=$(prefix)/lib/bem
foo_datafiles=\
	ic0.tri ic2.tri ic4.tri ic6.tri ic1.tri ic3.tri ic5.tri ic7.tri \
	inner_skull.dat outer_skull.dat outer_skin.dat

install-data-local:
	test -z $(DESTDIR)$(foodir) || $(mkdir_p) $(DESTDIR)$(foodir)
		$(INSTALL_DATA) $(foo_datafiles) $(DESTDIR)$(foodir)

# Our release target. Include files to be excluded here. They will be
# found and removed after 'make install' is run during the 'make
# release' target.
EXCLUDE_FILES=
include $(top_srcdir)/Makefile.extra
