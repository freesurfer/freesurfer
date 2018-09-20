project(vtkutils)

if(VTK_FOUND AND KWWidgets_FOUND)

  include_directories(
    ${FS_INCLUDE_DIRS}
    ${CMAKE_CURRENT_SOURCE_DIR}
    ${TIFF_INCLUDE_DIR}
    ${VTK_INCLUDE_DIRS}
    ${X11_INCLUDE_DIR}
    ${KWWidgets_INCLUDE_DIR}
  )

  set(CMAKE_CXX_STANDARD 98)

  # use vtkWrapTcl to create the tcl-wrapped code for these objects
  vtk_wrap_tcl(vtkKWRGBATransferFunctionEditor.h vtkKWRGBATransferFunctionEditorTcl.cxx)
  vtk_wrap_tcl(vtkKWBltGraph.h vtkKWBltGraphTcl.cxx)

  set(SOURCES
    vtkKWProgressDialog.cxx
    vtkKWRGBATransferFunctionEditor.cxx
    vtkKWRGBATransferFunctionEditorTcl.cxx
    vtkKWBltGraph.cxx vtkKWBltGraphTcl.cxx
    vtkFreesurferLookupTable.cxx
    vtkArrowPipeline.cxx
    vtkRGBATransferFunction.cxx
    vtkRGBAColorTransferFunction.cxx
    vtkFDTensorGlyph.cxx
    vtkODFGlyph.cxx
    vtkInflatePolyData.cxx
    IconLoader.cxx
  )

  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-deprecated")

  add_library(vtkutils STATIC ${SOURCES})
  target_link_libraries(vtkutils ${TIFF_LIBRARIES})

endif()
