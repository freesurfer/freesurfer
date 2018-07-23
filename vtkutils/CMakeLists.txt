project(vtkutils)

if(VTK_FOUND AND KWWidgets_FOUND AND TCLTKTIXBLT_FOUND)

  include_directories(
    ${CMAKE_SOURCE_DIR}/include
    ${TIFF_INCLUDE_DIR}
    ${VTK_INCLUDE_DIR}
    ${KWWidgets_INCLUDE_DIR}
  )

  # use vtkWrapTcl to create the tcl-wrapped code for these objects
  add_custom_command(
    OUTPUT vtkKWRGBATransferFunctionEditorTcl.cxx
    DEPENDS vtkKWRGBATransferFunctionEditor.h
    COMMAND ${VTK_WRAP_TCL_EXE} vtkKWRGBATransferFunctionEditor.h 
            ${VTK_LIBRARY_DIR}/hints 1 vtkKWRGBATransferFunctionEditorTcl.cxx
  )

  add_custom_command(
    OUTPUT vtkKWBltGraphTcl.cxx
    DEPENDS vtkKWBltGraph.h
    COMMAND ${VTK_WRAP_TCL_EXE} vtkKWBltGraph.h 
            ${VTK_LIBRARY_DIR}/hints 1 vtkKWBltGraphTcl.cxx
  )

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
  target_link_libraries(vtkutils ${TIFF_LIBRARIES} ${VTK_LIBRARIES} ${KWWidgets_LIBRARIES} ${TCLTKTIXBLT_LIBRARIES})

endif()
