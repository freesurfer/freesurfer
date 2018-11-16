project(vtkutils)

if(VTK_FOUND)

  include_directories(
    ${FS_INCLUDE_DIRS}
    ${CMAKE_CURRENT_SOURCE_DIR}
    ${TIFF_INCLUDE_DIR}
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

  # once qdec is deprecated, this can be removed entirely
  if(KWWidgets_FOUND AND BUILD_QDEC)
    include_directories(${KWWidgets_INCLUDE_DIR})

    # use vtkWrapTcl to create the tcl-wrapped code for these objects
    vtk_wrap_tcl(vtkKWRGBATransferFunctionEditor.h vtkKWRGBATransferFunctionEditorTcl.cxx)
    vtk_wrap_tcl(vtkKWBltGraph.h vtkKWBltGraphTcl.cxx)

    set(SOURCE ${SOURCES}
      vtkKWProgressDialog.cxx
      vtkKWRGBATransferFunctionEditor.cxx
      vtkKWRGBATransferFunctionEditorTcl.cxx
      vtkKWBltGraph.cxx
      vtkKWBltGraphTcl.cxx
      IconLoader.cxx
    )
  endif()

  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-deprecated")

  add_library(vtkutils STATIC ${SOURCES})
  target_link_libraries(vtkutils ${TIFF_LIBRARIES})

endif()
