project(histo_register)

if(OpenCV_FOUND)

  add_definitions(-DUSE_OPENCV -DUSE_CDT -DUSE_ZLIB)

  include_directories(
    ${FS_INCLUDE_DIRS}
    ${CMAKE_CURRENT_SOURCE_DIR}
    ${CMAKE_CURRENT_SOURCE_DIR}/SimpleBaseLib/include
    ${CMAKE_CURRENT_SOURCE_DIR}/SimpleBaseLib/external
    ${CMAKE_CURRENT_SOURCE_DIR}/ParticleVideoLib/include
    ${OpenCV_INCLUDE_DIR}
  )

  # local libraries
  add_subdirectory(SimpleBaseLib)
  add_subdirectory(ParticleVideoLib)

  set(SOURCES
    prep/BlockPrep.cc
    prep/HistoPrep.cc
    prep/HistoStats.cc
    prep/HistoStitch.cc
    prep/MPrep.cc
    prep/Polarization.cc
    prep/Prediction.cc
    prep/StitchNode.cc
    prep/VolumeFile.cc
    prep/VolumeUtil.cc
    registration/BlockRegister.cc
    registration/CorresField3D.cc
    registration/HistoRegister.cc
    registration/HistoTransform.cc
    registration/ImageSetSeq.cc
    registration/Normalization.cc
    registration/TestCorres3D.cc
    registration/VarCorres3D.cc
    registration/VarCorres3DUtil.cc
    Main.cc
  )

  add_executable(histo_register.bin ${SOURCES})
  target_link_libraries(histo_register.bin ${OpenCV_LIBRARIES} sbl pvl ${ZLIB_LIBRARIES})

  install(PROGRAMS histo_register_wrapper DESTINATION bin RENAME histo_register)
  install(TARGETS histo_register.bin DESTINATION bin)

endif()
