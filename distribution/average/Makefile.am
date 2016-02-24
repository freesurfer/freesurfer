#
# Makefile.am
#
SUBDIRS=surf \
	Yeo_JNeurophysiol11_MNI152 \
	Buckner_JNeurophysiol11_MNI152 \
	Choi_JNeurophysiol12_MNI152

foodir=$(prefix)/average

GCA=\
	face.gca \
	RB_all_2016-02-11.gca \
	RB_all_2015-08-04.gca \
	RB_all_2014-08-21.gca \
	RB_all_withskull_2014-08-21.gca \
	talairach_mixed_with_skull.gca \
	wmsa_new_eesmith.gca \
	aseg+spmhead.ixi.gca \
	aseg+spmhead+vermis+pons.ixi.gca \
	pons.mni152.2mm.mgz

GCS=\
	lh.atlas2002_simple.gcs \
	lh.atlas2005_simple.gcs \
	lh.curvature.buckner40.filled.desikan_killiany.2010-03-25.gcs \
	lh_trans_toSulc.gcs \
	lh.destrieux.simple.2009-07-29.gcs \
	rh.atlas2002_simple.gcs \
	rh.atlas2005_simple.gcs \
	rh.curvature.buckner40.filled.desikan_killiany.2010-03-25.gcs \
	rh_trans_toSulc.gcs \
	rh.destrieux.simple.2009-07-29.gcs \
	rh.DKTatlas40.gcs \
	lh.DKTatlas40.gcs \
	rh.DKTatlas100.gcs \
	lh.DKTatlas100.gcs

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
	colortable_BA.txt \
	tissue_parms.txt

MGZ=\
	HPMarch022005_ras.mgz \
	mni305.cor.mgz \
	mni305.cor.subfov1.mgz \
	mni305.cor.subfov2.mgz \
	mni305.mask.cor.mgz

if NMR_INSTALL
foo_datafiles=\
	$(GCA) \
	$(GCS) \
	$(TIF) \
	$(TXT) \
	$(MGZ) \
	label_scales.dat \
	mni305.cor.readme \
	mni305.cor.subfov1.reg \
	mni305.cor.subfov2.reg \
	mni152.register.dat \
	mni152.1mm.register.dat \
	mni152.mni305.cor.subfov1.dat \
	mni152.mni305.cor.subfov2.dat \
	RB_all_2008-03-26.mni152.2mm.lta 
else
foo_datafiles=\
	$(GCA) \
	$(GCS) \
	$(TIF) \
	$(TXT) \
	label_scales.dat \
	mni305.cor.mgz \
	mni305.cor.readme \
	mni305.mask.cor.mgz \
	mni305.cor.subfov1.mgz \
	mni305.cor.subfov2.mgz \
	mni305.cor.subfov1.reg \
	mni305.cor.subfov2.reg \
	mni152.register.dat \
	mni152.mni305.cor.subfov1.dat \
	mni152.mni305.cor.subfov2.dat \
	RB_all_2008-03-26.mni152.2mm.lta 
endif

foo_otherdatafiles=fsaverage.tar.gz fsaverage_sym.tar.gz \
        V1_average.tar.gz \
	fsaverage3.tar.gz fsaverage4.tar.gz \
	fsaverage5.tar.gz fsaverage6.tar.gz fsaverage7.tar.gz \
	EC_average.tar.gz mult-comp-cor.tar.gz Yeo_Brainmap_MNI152.tar.gz \
	build_gcs_morph.csh cvs_avg35.tar.gz cvs_avg35_inMNI152.tar.gz


install:
	test -z $(DESTDIR)$(foodir) || $(mkdir_p) $(DESTDIR)$(foodir)
	$(INSTALL_DATA) $(foo_datafiles) $(DESTDIR)$(foodir)
	mkdir -p $(DESTDIR)$(prefix)/subjects
	rm -Rf $(DESTDIR)$(prefix)/subjects/cvs_avg35
	rm -Rf $(top_builddir)/distribution/average/cvs_avg35
	gunzip -c $(top_srcdir)/distribution/average/cvs_avg35.tar.gz | tar xvf -
	mv $(top_builddir)/distribution/average/cvs_avg35 \
		$(DESTDIR)$(prefix)/subjects
	rm -Rf $(DESTDIR)$(prefix)/subjects/cvs_avg35_inMNI152
	rm -Rf $(top_builddir)/distribution/average/cvs_avg35_inMNI152
	gunzip -c $(top_srcdir)/distribution/average/cvs_avg35_inMNI152.tar.gz | tar xvf -
	mv $(top_builddir)/distribution/average/cvs_avg35_inMNI152 \
		$(DESTDIR)$(prefix)/subjects
	rm -Rf $(DESTDIR)$(prefix)/subjects/fsaverage
	rm -Rf $(top_builddir)/distribution/average/fsaverage
	gunzip -c $(top_srcdir)/distribution/average/fsaverage.tar.gz | tar xvf -
	mv $(top_builddir)/distribution/average/fsaverage \
		$(DESTDIR)$(prefix)/subjects
	chmod -R 775 $(DESTDIR)$(prefix)/subjects/fsaverage
	rm -Rf $(DESTDIR)$(prefix)/subjects/fsaverage_sym
	rm -Rf $(top_builddir)/distribution/average/fsaverage_sym
	gunzip -c $(top_srcdir)/distribution/average/fsaverage_sym.tar.gz | tar xvf -
	mv $(top_builddir)/distribution/average/fsaverage_sym \
		$(DESTDIR)$(prefix)/subjects
	chmod -R 775 $(DESTDIR)$(prefix)/subjects/fsaverage_sym
	rm -Rf $(DESTDIR)$(prefix)/subjects/V1_average
	rm -Rf $(top_builddir)/distribution/average/V1_average
	gunzip -c $(top_srcdir)/distribution/average/V1_average.tar.gz | tar xvf -
	cp -f $(top_builddir)/distribution/average/V1_average/scripts/predict_v1.sh \
		$(DESTDIR)$(bindir)
	mv $(top_builddir)/distribution/average/V1_average \
		$(DESTDIR)$(prefix)/subjects
	chmod -R 775 $(DESTDIR)$(prefix)/subjects/V1_average
	rm -Rf $(DESTDIR)$(prefix)/subjects/fsaverage3
	rm -Rf $(top_builddir)/distribution/average/fsaverage3
	gunzip -c $(top_srcdir)/distribution/average/fsaverage3.tar.gz | tar xvf -
	mv $(top_builddir)/distribution/average/fsaverage3 \
		$(DESTDIR)$(prefix)/subjects
	chmod -R 775 $(DESTDIR)$(prefix)/subjects/fsaverage3
	rm -Rf $(DESTDIR)$(prefix)/subjects/fsaverage4
	rm -Rf $(top_builddir)/distribution/average/fsaverage4
	gunzip -c $(top_srcdir)/distribution/average/fsaverage4.tar.gz | tar xvf -
	mv $(top_builddir)/distribution/average/fsaverage4 \
		$(DESTDIR)$(prefix)/subjects
	chmod -R 775 $(DESTDIR)$(prefix)/subjects/fsaverage4
	rm -Rf $(DESTDIR)$(prefix)/subjects/fsaverage5
	rm -Rf $(top_builddir)/distribution/average/fsaverage5
	gunzip -c $(top_srcdir)/distribution/average/fsaverage5.tar.gz | tar xvf -
	mv $(top_builddir)/distribution/average/fsaverage5 \
		$(DESTDIR)$(prefix)/subjects
	chmod -R 775 $(DESTDIR)$(prefix)/subjects/fsaverage5
	rm -Rf $(DESTDIR)$(prefix)/subjects/fsaverage6
	rm -Rf $(top_builddir)/distribution/average/fsaverage6
	gunzip -c $(top_srcdir)/distribution/average/fsaverage6.tar.gz | tar xvf -
	mv $(top_builddir)/distribution/average/fsaverage6 \
		$(DESTDIR)$(prefix)/subjects
	chmod -R 775 $(DESTDIR)$(prefix)/subjects/fsaverage6
	rm -Rf $(DESTDIR)$(prefix)/subjects/?h.EC_average
	rm -Rf $(top_builddir)/distribution/average/?h.EC_average
	gunzip -c $(top_srcdir)/distribution/average/EC_average.tar.gz | tar xvf -
	mv $(top_builddir)/distribution/average/?h.EC_average \
		$(DESTDIR)$(prefix)/subjects
	chmod -R 775 $(DESTDIR)$(prefix)/subjects/?h.EC_average
	rm -Rf $(DESTDIR)$(prefix)/average/mult-comp-cor
	rm -Rf $(top_builddir)/distribution/average/mult-comp-cor
	gunzip -c $(top_srcdir)/distribution/average/mult-comp-cor.tar.gz | tar xvf -
	pwd
	mv $(top_builddir)/distribution/average/mult-comp-cor \
		$(DESTDIR)$(prefix)/average
	chmod -R 775 $(DESTDIR)$(prefix)/average/mult-comp-cor
	rm -Rf $(DESTDIR)$(prefix)/average/Yeo_Brainmap_MNI152
	rm -Rf $(top_builddir)/distribution/average/Yeo_Brainmap_MNI152
	gunzip -c $(top_srcdir)/distribution/average/Yeo_Brainmap_MNI152.tar.gz | tar xvf -
	pwd
	mv $(top_builddir)/distribution/average/Yeo_Brainmap_MNI152 \
		$(DESTDIR)$(prefix)/average
	chmod -R 775 $(DESTDIR)$(prefix)/average/Yeo_Brainmap_MNI152
	cd surf; make install
	cd Yeo_JNeurophysiol11_MNI152; make install
	cd Buckner_JNeurophysiol11_MNI152; make install
	cd Choi_JNeurophysiol12_MNI152; make install

uninstall-hook:
	rm -Rf $(DESTDIR)$(prefix)/subjects/cvs_avg35
	rm -Rf $(DESTDIR)$(prefix)/subjects/cvs_avg35_inMNI152
	rm -Rf $(DESTDIR)$(prefix)/subjects/fsaverage
	rm -Rf $(DESTDIR)$(prefix)/subjects/fsaverage_sym
	rm -Rf $(DESTDIR)$(prefix)/subjects/fsaverage3
	rm -Rf $(DESTDIR)$(prefix)/subjects/fsaverage4
	rm -Rf $(DESTDIR)$(prefix)/subjects/fsaverage5
	rm -Rf $(DESTDIR)$(prefix)/subjects/fsaverage6
	rm -Rf $(DESTDIR)$(prefix)/subjects/?h.EC_average
	rm -Rf $(DESTDIR)$(prefix)/subjects/V1_average
	rm -f $(DESTDIR)$(bindir)/predict_v1.sh
	rm -Rf $(DESTDIR)$(prefix)/average/mult-comp-cor
	rm -Rf $(DESTDIR)$(prefix)/average/Yeo_Brainmap_MNI152

# Our release target. Include files to be excluded here. They will be
# found and removed after 'make install' is run during the 'make
# release' target.
EXCLUDE_FILES="HPMarch022005_ras.mgz"

include $(top_srcdir)/Makefile.extra
