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
 $(wildcard *.mgz) \
 mni305.cor.readme 
else
foo_DATA=$(wildcard *.gcs) \
 $(wildcard *.gca) \
 $(wildcard *.tif) \
 $(wildcard *.txt) \
 mni305.cor.mgz \
 mni305.cor.readme \
 mni305.mask.cor.mgz
endif

EXTRA_DIST=$(foo_DATA)

install: install-fooDATA
	mkdir -p $(prefix)/subjects
	tar zxvf fsaverage.tar.gz --modification-time --directory $(prefix)/subjects

# Our release target. Include files to be excluded here. They will be
# found and removed after 'make install' is run during the 'make
# release' target.
EXCLUDE_FILES="HPMarch022005_ras.mgz"
include $(top_srcdir)/Makefile.extra

# for the make dist target, which creates the source-code tarball,
# use 'dist-hook' to do what EXCLUDE_FILES does for the make install target,
# that is, remove files that should not go public.  
dist-hook:
	rm -f $(distdir)/distribution/average/HPMarch022005_ras.mgz
