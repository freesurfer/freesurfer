#
# Makefile.am
#
SUBDIRS=surf

foodir=$(prefix)/average

GCA=\
	face.gca \
	RB_all_2008-03-26.gca \
	RB_all_withskull_2008-03-26.gca

GCS=\
	lh.atlas2002_simple.gcs \
	lh.atlas2005_simple.gcs \
	lh.curvature.buckner40.filled.desikan_killiany.2007-06-20.gcs \
	lh_trans_toSulc.gcs \
	rh.atlas2002_simple.gcs \
	rh.atlas2005_simple.gcs \
	rh.curvature.buckner40.filled.desikan_killiany.2007-06-20.gcs \
	rh_trans_toSulc.gcs

TIF=\
	lh.average.curvature.filled.buckner40.tif \
	lh.average.CURVATURE.tif \
	lh.average.tif \
	rh.average.curvature.filled.buckner40.tif \
	rh.average.CURVATURE.tif \
	rh.average.tif \
	rigidly_aligned_brain_template.tif

TXT=\
	colortable_desikan_killiany.txt \
	tissue_parms.txt

MGZ=\
	HPMarch022005_ras.mgz \
	mni305.cor.mgz \
	mni305.cor.subfov1.mgz \
	mni305.cor.subfov2.mgz \
	mni305.mask.cor.mgz

if NMR_INSTALL
foo_DATA=\
	$(GCA) \
	$(GCS) \
	$(TIF) \
	$(TXT) \
	$(MGZ) \
	mni305.cor.readme \
	mni305.cor.subfov1.reg \
	mni305.cor.subfov2.reg \
	mni152.register.dat
else
foo_DATA=\
	$(GCA) \
	$(GCS) \
	$(TIF) \
	$(TXT) \
	mni305.cor.mgz \
	mni305.cor.readme \
	mni305.mask.cor.mgz \
	mni305.cor.subfov1.mgz \
	mni305.cor.subfov2.mgz \
	mni305.cor.subfov1.reg \
	mni305.cor.subfov2.reg \
	mni152.register.dat
endif

EXTRA_DIST=$(foo_DATA) fsaverage.tar.gz V1_average.tar.gz

install: install-fooDATA
	mkdir -p $(DESTDIR)$(prefix)/subjects
	rm -Rf $(DESTDIR)$(prefix)/subjects/fsaverage
	gunzip -c $(top_srcdir)/distribution/average/fsaverage.tar.gz | tar xvf -
	mv $(top_builddir)/distribution/average/fsaverage \
		$(DESTDIR)$(prefix)/subjects
	rm -Rf $(DESTDIR)$(prefix)/subjects/V1_average
	gunzip -c $(top_srcdir)/distribution/average/V1_average.tar.gz | tar xvf -
	mv $(top_builddir)/distribution/average/V1_average \
		$(DESTDIR)$(prefix)/subjects
	cd surf; make install

uninstall-hook:
	rm -Rf $(DESTDIR)$(prefix)/subjects/fsaverage
	rm -Rf $(DESTDIR)$(prefix)/subjects/V1_average

# Our release target. Include files to be excluded here. They will be
# found and removed after 'make install' is run during the 'make
# release' target.
EXCLUDE_FILES="HPMarch022005_ras.mgz"

include $(top_srcdir)/Makefile.extra
