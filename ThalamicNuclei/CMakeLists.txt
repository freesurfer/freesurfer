project(ThalamicNuclei)

if(APPLE)
  set(BIN_PROGRAMS
    mac_osx/run_SegmentThalamicNuclei.sh
  )
  # todo: install the whole app!
else()
  set(BIN_PROGRAMS
    linux_x86_64/SegmentThalamicNuclei
    linux_x86_64/run_SegmentThalamicNuclei.sh
  )
endif()

install_symlinks(
  ${BIN_PROGRAMS}
  segmentThalamicNuclei.sh
  quantifyThalamicNuclei.sh
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