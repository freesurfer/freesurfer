# these are the minimum-required files for a standard recon
install_symlinks(
  TYPE files
  DESTINATION average
  # gcs
  lh.DKaparc.atlas.acfb40.noaparc.i12.2016-08-02.gcs
  rh.DKaparc.atlas.acfb40.noaparc.i12.2016-08-02.gcs
  lh.DKTaparc.atlas.acfb40.noaparc.i12.2020-05-13.gcs
  rh.DKTaparc.atlas.acfb40.noaparc.i12.2020-05-13.gcs
  lh.DKTaparc.atlas.acfb40.noaparc.i12.2016-08-02.gcs
  rh.DKTaparc.atlas.acfb40.noaparc.i12.2016-08-02.gcs
  lh.CDaparc.atlas.acfb40.noaparc.i12.2016-08-02.gcs
  rh.CDaparc.atlas.acfb40.noaparc.i12.2016-08-02.gcs
  # gca
  RB_all_2016-05-10.vc700.gca
  RB_all_withskull_2016-05-10.vc700.gca
  RB_all_2019_10_25.talxfm.mni305.gca
  RB_all_withskull_2019_10_22.talxfm.mni305.gca
  RB_all_2020-01-02.gca
  RB_all_withskull_2020_01_02.gca
  # tiff
  rigidly_aligned_brain_template.tif
  lh.folding.atlas.acfb40.noaparc.i12.2016-08-02.tif
  rh.folding.atlas.acfb40.noaparc.i12.2016-08-02.tif
  # txt
  colortable_desikan_killiany.txt
  colortable_BA.txt
  colortable_vpnl.txt
  # mgz
  mni305.cor.mgz
  mni305.cor.subfov1.mgz
  mni305.cor.subfov2.mgz
  mni305.mask.cor.mgz
)

# fsaverage
install_tarball(fsaverage.tar.gz DESTINATION subjects)

# samseg atlas
add_subdirectory(samseg)

# subjects required for infant recon-all
if(INFANT_MODULE)
  install_tarball(CNYBCH.tar.gz DESTINATION average)
endif()

if(INFANT_MODULE AND MINIMAL)
  install_symlinks(
    TYPE files
    DESTINATION average
    lh.average.curvature.filled.buckner40.tif
    rh.average.curvature.filled.buckner40.tif
    lh.curvature.buckner40.filled.desikan_killiany.2010-03-25.gcs
    rh.curvature.buckner40.filled.desikan_killiany.2010-03-25.gcs
    lh.destrieux.simple.2009-07-29.gcs
    rh.destrieux.simple.2009-07-29.gcs
  )
endif()

# the following data isn't required for recon-all, but will be installed by default
if(NOT MINIMAL)

  add_subdirectories(
    surf
    Yeo_JNeurophysiol11_MNI152
    Buckner_JNeurophysiol11_MNI152
    Choi_JNeurophysiol12_MNI152
    mideface-atlas
  )

  install_symlinks(
    TYPE files
    DESTINATION average
    # ca
    face.gca
    talairach_mixed_with_skull.gca
    wmsa_new_eesmith.gca
    aseg+spmhead+vermis+pons.ixi.gca
    pons.mni152.2mm.mgz
    # gcs
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
    # tif
    lh.average.curvature.filled.buckner40.tif
    lh.average.CURVATURE.tif
    lh.average.tif
    rh.average.curvature.filled.buckner40.tif
    rh.average.CURVATURE.tif
    rh.average.tif
    # dat
    label_scales.dat
    mni305.cor.readme
    mni305.cor.subfov1.reg
    mni305.cor.subfov2.reg
    mni152.register.dat
    mni152.mni305.cor.subfov1.dat
    mni152.mni305.cor.subfov2.dat
    RB_all_2008-03-26.mni152.2mm.lta
    # txt
    tissue_parms.txt
  )

  # nmr-only files
  install_symlinks(
    HPMarch022005_ras.mgz
    aseg+spmhead.ixi.gca
    mni152.1mm.register.dat
    TYPE files
    DESTINATION average
    NMR_ONLY
  )

  # unpack subjects
  install(PROGRAMS predict_v1.sh DESTINATION bin)
  install_tarball(cvs_avg35.tar.gz DESTINATION subjects)
  install_tarball(cvs_avg35_inMNI152.tar.gz DESTINATION subjects)
  install_tarball(EC_average.tar.gz DESTINATION subjects)
  install_tarball(fsaverage3.tar.gz DESTINATION subjects)
  install_tarball(fsaverage4.tar.gz DESTINATION subjects)
  install_tarball(fsaverage5.tar.gz DESTINATION subjects)
  install_tarball(fsaverage6.tar.gz DESTINATION subjects)
  install_tarball(fsaverage_sym.tar.gz DESTINATION subjects)
  install_tarball(V1_average.tar.gz DESTINATION subjects)
  install_tarball(mult-comp-cor.tar.gz DESTINATION average)
  install_tarball(Yeo_Brainmap_MNI152.tar.gz DESTINATION average)

endif()
