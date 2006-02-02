#
# Makefile.am
#

SUBDIRS=surf

foodir=$(prefix)/average
if NMR_INSTALL
foo_DATA=$(wildcard *.gcs) \
 $(wildcard *.gca) \
 $(wildcard *.tif) \
 $(wildcard *.txt) \
 $(wildcard *.mgz)
else
foo_DATA=$(wildcard *.gcs) \
 $(wildcard *.gca) \
 $(wildcard *.tif) \
 $(wildcard *.txt)
endif

EXTRA_DIST=$(foo_DATA)

# Our release target. Include files to be excluded here. They will be
# found and removed after 'make install' is run during the 'make
# release' target.
EXCLUDE_FILES="HPMarch022005_ras.mgz"
include $(top_srcdir)/Makefile.extra

# for the make dist target, which creates the source-code tarball,
# use 'dist-hook' to do what EXCLUDE_FILES does for the make install target,
# that is, remove files that should not go public.  
dist-hook:
	rm -f HPMarch022005_ras.mgz
