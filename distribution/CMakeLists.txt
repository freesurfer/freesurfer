install(FILES
  ASegStatsLUT.txt
  DefectLUT.txt
  FreeSurferColorLUT.txt
  numaris4-protocols.unpackcfg
  SegmentNoLUT.txt
  Simple_surface_labels2005.txt
  Simple_surface_labels2009.txt
  tkmeditParcColorsCMA
  WMParcStatsLUT.txt
  DESTINATION ${CMAKE_INSTALL_PREFIX}
)

install(PROGRAMS
  FreeSurferEnv.csh
  FreeSurferEnv.sh
  SetUpFreeSurfer.csh
  SetUpFreeSurfer.sh
  sources.csh
  sources.sh
  DESTINATION ${CMAKE_INSTALL_PREFIX}
)

set(FSLDIR $ENV{FSLDIR})
if(FSLDIR)
  foreach(CMD flirt fslmaths slicetimer bet fslswapdim fslorient)
    install(PROGRAM $(FSLDIR)/bin/${CMD} DESTINATION bin RENAME ${CMD}.fsl)
  endforeach()
endif()

set(DCMTKDIR $ENV{DCMTKDIR})
if(DCMTKDIR)
  install(PROGRAM $(FSLDIR)/bin/dcmdjpeg DESTINATION bin RENAME dcmdjpeg.fs)
  install(PROGRAM $(FSLDIR)/bin/dcmdrle DESTINATION bin RENAME dcmdrle.fs)
endif()

add_subdirectory(average)
add_subdirectory(diffusion)
add_subdirectory(docs)
add_subdirectory(fsafd)
add_subdirectory(subjects)

install(FILES sessions/README DESTINATION sessions)
install_tarball(trctrain/trctraindata.tar.gz trctrain)
