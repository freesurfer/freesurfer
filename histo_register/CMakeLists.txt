project(histo_register)

if(OpenCV_FOUND)

  add_definitions(-DUSE_OPENCV -DUSE_CDT -DUSE_ZLIB)

  include_directories(
    ${CMAKE_SOURCE_DIR}/include
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

  set(WRAPCODE "\
#!/bin/tcsh -ef\n\
setenv LD_LIBRARY_PATH /usr/pubsw/packages/opencv/current/lib\n\
histo_register.bin $argv"
  )
  add_executable(histo_register ${SOURCES})
  target_link_libraries(histo_register ${OpenCV_LIBRARIES} sbl pvl)
  install_wrapped(TARGETS histo_register DESTINATION bin WRAPCODE ${WRAPCODE})

endif()
