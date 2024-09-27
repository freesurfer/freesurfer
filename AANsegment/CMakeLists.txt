project(AANsegment)

if(APPLE)
  set(BIN_PROGRAMS
    mac_osx/run_segmentNuclei.sh
  )
  install_osx_app(mac_osx/segmentNuclei.app)
else()
  set(BIN_PROGRAMS
    linux_x86_64/segmentNuclei
    linux_x86_64/run_segmentNuclei.sh
  )
endif()

install_symlinks(
  ${BIN_PROGRAMS}
  SegmentAAN.sh
  TYPE programs
  DESTINATION bin
)

install_symlinks(
  targetReg.mgz
  targetWorkingres.mgz
  AtlasMesh.gz
  compressionLookupTable.txt
  freeview.lut.txt
  TYPE files
  DESTINATION average/AAN/atlas
)
