#
# Makefile.am
#

SUBDIRS=surf

foodir=$(prefix)/average
foo_DATA=$(wildcard *.gcs) \
 $(wildcard *.gca) \
 $(wildcard *.tif) \
 $(wildcard *.txt) \
 $(wildcard *.ctab) \
 $(wildcard *.mgz)

EXTRA_DIST=$(foo_DATA)

# Our release target. Include files to be excluded here. They will be
# found and removed after 'make install' is run during the 'make
# release' target.
EXCLUDE_FILES="HPMarch022005_ras.mgz"
include $(top_srcdir)/Makefile.extra
