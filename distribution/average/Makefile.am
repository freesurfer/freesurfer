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
 mni305.cor.readme \
 mni305.cor.subfov1.reg \
 mni305.cor.subfov2.reg \
 mni152.register.dat 
else
foo_DATA=$(wildcard *.gcs) \
 $(wildcard *.gca) \
 $(wildcard *.tif) \
 $(wildcard *.txt) \
 mni305.cor.mgz \
 mni305.cor.readme \
 mni305.mask.cor.mgz \
 mni305.cor.subfov1.mgz \
 mni305.cor.subfov2.mgz \
 mni305.cor.subfov1.reg \
 mni305.cor.subfov2.reg \
 mni152.register.dat 
endif

EXTRA_DIST=$(foo_DATA) fsaverage.tar.gz

install: install-fooDATA
	mkdir -p $(DESTDIR)$(prefix)/subjects
	rm -Rf $(DESTDIR)$(prefix)/subjects/fsaverage
	gunzip -c fsaverage.tar.gz | tar xvf -
	mv fsaverage $(DESTDIR)$(prefix)/subjects

# Our release target. Include files to be excluded here. They will be
# found and removed after 'make install' is run during the 'make
# release' target.
EXCLUDE_FILES="HPMarch022005_ras.mgz"
include $(top_srcdir)/Makefile.extra
