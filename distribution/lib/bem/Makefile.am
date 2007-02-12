#
# Makefile.am
#
foodir=$(prefix)/lib/bem
foo_DATA=$(wildcard *.tri) $(wildcard *.dat)

EXTRA_DIST=$(foo_DATA)

# Our release target. Include files to be excluded here. They will be
# found and removed after 'make install' is run during the 'make
# release' target.
EXCLUDE_FILES=
include $(top_srcdir)/Makefile.extra
