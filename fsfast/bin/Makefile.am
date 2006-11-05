#
# Makefile.am for fsfast/bin
#

foodir=$(prefix)/fsfast/bin
foo_SCRIPTS=$(wildcard *)

# foo_SCRIPTS has all-files wildcard, but lets not install Makefile stuff.
# note: install-exec-hook is already defined in Makefile.extra
install-data-hook:
	rm -f $(DESTDIR)$(foodir)/Makefile*

EXTRA_DIST=$(foo_SCRIPTS)

# Our release target. Include files to be excluded here. They will be
# found and removed after 'make install' is run during the 'make
# release' target.
TDR_FILES=$(wildcard tdr-*)
EXCLUDE_FILES=$(TDR_FILES) CVS
include $(top_srcdir)/Makefile.extra
