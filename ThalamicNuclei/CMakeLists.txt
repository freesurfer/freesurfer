project(ThalamicNuclei)

if(APPLE)
  set(BIN_PROGRAMS
    mac_osx/run_SegmentThalamicNuclei.sh
  )
  install_osx_app(mac_osx/SegmentThalamicNuclei.app)
else()
  set(BIN_PROGRAMS
    linux_x86_64/SegmentThalamicNuclei
    linux_x86_64/run_SegmentThalamicNuclei.sh
  )
endif()

install_symlinks(
  ${BIN_PROGRAMS}
  segmentThalamicNuclei.sh
  TYPE programs
  DESTINATION bin
)

install_symlinks(
  AtlasDump.mgz
  AtlasMesh.gz
  compressionLookupTable.txt
  TYPE files
  DESTINATION average/ThalamicNuclei/atlas
)