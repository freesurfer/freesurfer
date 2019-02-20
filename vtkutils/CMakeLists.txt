project(vtkutils)

if(VTK_FOUND)

  include_directories(
    ${FS_INCLUDE_DIRS}
    ${CMAKE_CURRENT_SOURCE_DIR}
    ${CMAKE_SOURCE_DIR}/packages/tiff
    ${VTK_INCLUDE_DIRS}
    ${X11_INCLUDE_DIR}
  )

  set(CMAKE_CXX_STANDARD 98)

  set(SOURCES
    vtkFreesurferLookupTable.cxx
    vtkArrowPipeline.cxx
    vtkRGBATransferFunction.cxx
    vtkRGBAColorTransferFunction.cxx
    vtkFDTensorGlyph.cxx
    vtkODFGlyph.cxx
    vtkInflatePolyData.cxx
  )

  if(SUPPRESS_WARNINGS)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-deprecated")
  endif()

  add_library(vtkutils STATIC ${SOURCES})
  target_link_libraries(vtkutils tiff)

endif()
