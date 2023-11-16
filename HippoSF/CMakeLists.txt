project(HippoSF)

if(APPLE)
  set(BIN_PROGRAMS
    mac_osx/run_segmentSubjectT1_autoEstimateAlveusML.sh
    mac_osx/run_segmentSubjectT1T2_autoEstimateAlveusML.sh
    mac_osx/run_segmentSubjectT2_autoEstimateAlveusML.sh
    mac_osx/run_SegmentSubfieldsT1Longitudinal.sh
  )
  install_osx_app(mac_osx/SegmentSubfieldsT1Longitudinal.app)
  install_osx_app(mac_osx/segmentSubjectT1_autoEstimateAlveusML.app)
  install_osx_app(mac_osx/segmentSubjectT1T2_autoEstimateAlveusML.app)
  install_osx_app(mac_osx/segmentSubjectT2_autoEstimateAlveusML.app)
else()
  set(BIN_PROGRAMS
    linux_x86_64/run_SegmentSubfieldsT1Longitudinal.sh
    linux_x86_64/run_segmentSubjectT1_autoEstimateAlveusML.sh
    linux_x86_64/run_segmentSubjectT1T2_autoEstimateAlveusML.sh
    linux_x86_64/run_segmentSubjectT2_autoEstimateAlveusML.sh
    linux_x86_64/SegmentSubfieldsT1Longitudinal
    linux_x86_64/segmentSubjectT1_autoEstimateAlveusML
    linux_x86_64/segmentSubjectT1T2_autoEstimateAlveusML
    linux_x86_64/segmentSubjectT2_autoEstimateAlveusML
  )
endif()

install_symlinks(
  ${BIN_PROGRAMS}
  segmentHA_T1_long.sh
  segmentHA_T1.sh
  segmentHA_T2.sh
  TYPE programs
  DESTINATION bin
)

install_symlinks(
  AtlasDump.mgz
  AtlasMesh.gz
  compressionLookupTable.txt
  TYPE files
  DESTINATION average/HippoSF/atlas
)
