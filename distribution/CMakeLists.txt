install(FILES
  ASegStatsLUT.txt
  DefectLUT.txt
  FreeSurferColorLUT.txt
  SegmentNoLUT.txt
  Simple_surface_labels2005.txt
  Simple_surface_labels2009.txt
  SubCorticalMassLUT.txt
  tkmeditParcColorsCMA
  WMParcStatsLUT.txt
  DESTINATION .
)

install_configured(
  FreeSurferEnv.csh
  FreeSurferEnv.sh
  DESTINATION .
)

install(PROGRAMS
  SetUpFreeSurfer.csh
  SetUpFreeSurfer.sh
  sources.csh
  sources.sh
  DESTINATION .
)

add_subdirectories(average docs fsafd subjects etc)

install_symlinks(
  lib/bem/ic0.tri
  lib/bem/ic1.tri
  lib/bem/ic2.tri
  lib/bem/ic3.tri
  lib/bem/ic4.tri
  lib/bem/ic5.tri
  lib/bem/ic6.tri
  lib/bem/ic7.tri
  lib/bem/inner_skull.dat
  lib/bem/outer_skin.dat
  lib/bem/outer_skull.dat
  TYPE files
  DESTINATION lib/bem
)

# these aren't required for recon-all, but will be installed by default
if(NOT MINIMAL)
  add_subdirectory(diffusion)
  install(FILES sessions/README DESTINATION sessions)
  install_tarball(trctrain.tar.gz DESTINATION trctrain)
endif()

# FSL commands
set(FSLDIR $ENV{FSLDIR})
if(FSLDIR)
  foreach(CMD flirt fslmaths slicetimer bet fslswapdim fslorient)
    install(PROGRAMS ${FSLDIR}/bin/${CMD} DESTINATION bin RENAME ${CMD}.fsl)
  endforeach()
endif()

# DCMTK commands
set(DCMTKDIR $ENV{DCMTKDIR})
if(DCMTKDIR)
  install(PROGRAMS ${DCMTKDIR}/bin/dcmdjpeg DESTINATION bin RENAME dcmdjpeg.fs)
  install(PROGRAMS ${DCMTKDIR}/bin/dcmdrle DESTINATION bin RENAME dcmdrle.fs)
endif()

# GDCM commands
set(GDCMDIR $ENV{GDCMDIR})
if(GDCMDIR)
  install(PROGRAMS ${GDCMDIR}/bin/gdcmconv DESTINATION bin RENAME gdcmconv.fs)
endif()

# AFNI commands
set(AFNIDIR $ENV{AFNIDIR})
if(AFNIDIR)
  install(PROGRAMS ${AFNIDIR}/bin/3dvolreg DESTINATION bin RENAME 3dvolreg.afni)
endif()

# symlink MNI installation
if(EXISTS ${FS_PACKAGES_DIR}/mni/current)
  symlink(${FS_PACKAGES_DIR}/mni/current ${CMAKE_INSTALL_PREFIX}/mni)
endif()
