#
# Makefile.am for matlab files
#
# here even though file list is not given, but .m is enough to discriminate files
#
foodir=$(prefix)/matlab
foo_DATA= $(wildcard *.m)

# Our release target. Include files to be excluded here. They will be
# found and removed after 'make install' is run during the 'make
# release' target.
EXCLUDE_FILES=""
include $(top_srcdir)/Makefile.extra
