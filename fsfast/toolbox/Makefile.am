#
# Makefile.am for fsfast/toolbox
#

foodir=$(prefix)/fsfast/toolbox
foo_DATA= $(wildcard *.m *.fig)

EXTRA_DIST=$(foo_DATA)

# Our release target. Include files to be excluded here. They will be
# found and removed after 'make install' is run during the 'make
# release' target.
TDR_FILES=$(wildcard tdr_*.m)
EXCLUDE_FILES=$(TDR_FILES)
include $(top_srcdir)/Makefile.extra
