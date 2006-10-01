#
# Makefile.am for fsfast/bin
#

foodir=$(prefix)/fsfast/bin
foo_SCRIPTS=$(wildcard *)

EXTRA_DIST=$(foo_SCRIPTS)

# Our release target. Include files to be excluded here. They will be
# found and removed after 'make install' is run during the 'make
# release' target.
TDR_FILES=$(wildcard tdr-*)
EXCLUDE_FILES=$(TDR_FILES)
include $(top_srcdir)/Makefile.extra
