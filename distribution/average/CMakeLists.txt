# CA
install_symlinks(
  face.gca
  RB_all_2016-05-10.vc700.gca
  RB_all_withskull_2016-05-10.vc700.gca
  talairach_mixed_with_skull.gca
  wmsa_new_eesmith.gca
  aseg+spmhead.ixi.gca
  aseg+spmhead+vermis+pons.ixi.gca
  pons.mni152.2mm.mgz
  TYPE FILES
  DESTINATION average
)

# GCS
install_symlinks(
  lh.atlas2002_simple.gcs
  lh.atlas2005_simple.gcs
  lh.curvature.buckner40.filled.desikan_killiany.2010-03-25.gcs
  lh_trans_toSulc.gcs
  lh.destrieux.simple.2009-07-29.gcs
  rh.atlas2002_simple.gcs
  rh.atlas2005_simple.gcs
  rh.curvature.buckner40.filled.desikan_killiany.2010-03-25.gcs
  rh_trans_toSulc.gcs
  rh.destrieux.simple.2009-07-29.gcs
  rh.DKTatlas40.gcs
  lh.DKTatlas40.gcs
  rh.DKTatlas100.gcs
  lh.DKTatlas100.gcs
  lh.DKaparc.atlas.acfb40.noaparc.i12.2016-08-02.gcs
  rh.DKaparc.atlas.acfb40.noaparc.i12.2016-08-02.gcs
  lh.DKTaparc.atlas.acfb40.noaparc.i12.2016-08-02.gcs
  rh.DKTaparc.atlas.acfb40.noaparc.i12.2016-08-02.gcs
  lh.CDaparc.atlas.acfb40.noaparc.i12.2016-08-02.gcs
  rh.CDaparc.atlas.acfb40.noaparc.i12.2016-08-02.gcs
  TYPE FILES
  DESTINATION average
)

# TIF
install_symlinks(
  lh.average.curvature.filled.buckner40.tif
  lh.average.CURVATURE.tif
  lh.average.tif
  rh.average.curvature.filled.buckner40.tif
  rh.average.CURVATURE.tif
  rh.average.tif
  rigidly_aligned_brain_template.tif
  lh.folding.atlas.acfb40.noaparc.i12.2016-08-02.tif
  rh.folding.atlas.acfb40.noaparc.i12.2016-08-02.tif
  TYPE FILES
  DESTINATION average
)

# TXT
install_symlinks(
  colortable_desikan_killiany.txt
  colortable_BA.txt
  tissue_parms.txt
  TYPE FILES
  DESTINATION average
)

# MGZ
install_symlinks(
  mni305.cor.mgz
  mni305.cor.subfov1.mgz
  mni305.cor.subfov2.mgz
  mni305.mask.cor.mgz
  TYPE FILES
  DESTINATION average
)

# DAT
install_symlinks(
  label_scales.dat
  mni305.cor.readme
  mni305.cor.subfov1.reg
  mni305.cor.subfov2.reg
  mni152.register.dat
  mni152.mni305.cor.subfov1.dat
  mni152.mni305.cor.subfov2.dat
  RB_all_2008-03-26.mni152.2mm.lta
  TYPE FILES
  DESTINATION average
)

# unpack subjects
install_tarball(fsaverage.tar.gz subjects)
install_tarball(cvs_avg35.tar.gz subjects)
install_tarball(cvs_avg35_inMNI152.tar.gz subjects)
install_tarball(fsaverage_sym.tar.gz subjects)
install_tarball(V1_average.tar.gz subjects)
install_tarball(fsaverage3.tar.gz subjects)
install_tarball(fsaverage4.tar.gz subjects)
install_tarball(fsaverage5.tar.gz subjects)
install_tarball(fsaverage6.tar.gz subjects)
install_tarball(EC_average.tar.gz subjects)
install_tarball(mult-comp-cor.tar.gz average)
install_tarball(Yeo_Brainmap_MNI152.tar.gz average)

# subdirectories
add_subdirectory(surf)
add_subdirectory(samseg)
add_subdirectory(Yeo_JNeurophysiol11_MNI152)
add_subdirectory(Buckner_JNeurophysiol11_MNI152)
add_subdirectory(Choi_JNeurophysiol12_MNI152)
