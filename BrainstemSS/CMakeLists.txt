project(BrainstemSS)

if(APPLE)
  set(BIN_PROGRAMS mac_osx/run_SegmentSubject.sh)
  install_osx_app(mac_osx/segmentSubject.app)
else()
  set(BIN_PROGRAMS
    linux_x86_64/run_SegmentSubject.sh
    linux_x86_64/segmentSubject
  )
endif()

install_symlinks(
  ${BIN_PROGRAMS}
  segmentBS.sh
  ConcatenateSubregionsResults.sh
  TYPE programs
  DESTINATION bin
)

install_symlinks(
  AtlasDump.mgz
  AtlasMesh.gz
  compressionLookupTable.txt
  TYPE files
  DESTINATION average/BrainstemSS/atlas
)
