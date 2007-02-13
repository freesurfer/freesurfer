#
# Makefile.am
#

foodir=$(prefix)/average/surf
foo_DATA=lh.sphere.reg rh.sphere.reg

# do nothing when making distclean:
distclean:

EXTRA_DIST=$(foo_DATA)

# Our release target. Include files to be excluded here. They will be
# found and removed after 'make install' is run during the 'make
# release' target.
EXCLUDE_FILES=
include $(top_srcdir)/Makefile.extra
