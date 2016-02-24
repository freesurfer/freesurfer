#
# Makefile.am
#

foodir=$(prefix)/average/surf

foo_datafiles=README lh.sphere.reg rh.sphere.reg \
lh.sphere.ico0.reg  lh.sphere.ico4.reg  rh.sphere.ico0.reg  rh.sphere.ico4.reg \
lh.sphere.ico1.reg  lh.sphere.ico5.reg  rh.sphere.ico1.reg  rh.sphere.ico5.reg \
lh.sphere.ico2.reg  lh.sphere.ico6.reg  rh.sphere.ico2.reg  rh.sphere.ico6.reg \
lh.sphere.ico3.reg  lh.sphere.ico7.reg  rh.sphere.ico3.reg  rh.sphere.ico7.reg

install-data-local:
	test -z $(DESTDIR)$(foodir) || $(mkdir_p) $(DESTDIR)$(foodir)
	$(INSTALL_DATA) $(foo_datafiles) $(DESTDIR)$(foodir)

# Our release target. Include files to be excluded here. They will be
# found and removed after 'make install' is run during the 'make
# release' target.
EXCLUDE_FILES=
include $(top_srcdir)/Makefile.extra
