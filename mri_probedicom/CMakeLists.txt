project(mri_probedicom)

if(OPENGL_FOUND)

  include_directories(
    ${FS_INCLUDE_DIRS}
    ${CMAKE_SOURCE_DIR}/packages/dicom
    ${CMAKE_SOURCE_DIR}/packages/glut
    ${CMAKE_SOURCE_DIR}/packages/dcm2niix
    ${X11_INCLUDE_DIR}
  )

  add_executable(mri_probedicom mri_probedicom.cpp)
  if(${CMAKE_CXX_COMPILER_ID} STREQUAL "Clang")
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -Wl,-dead_strip -Wl,-stack_size -Wl,0x1000000")
  endif()
  target_link_libraries(mri_probedicom utils ${OPENGL_LIBRARIES} glut ${X11_LIBRARIES})
  install(TARGETS mri_probedicom DESTINATION bin)

endif()
