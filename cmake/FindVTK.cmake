# VTK Find Module

# VTK6+ has a useful cmake config that sets up
# variables and component paths nicely. Unfortunately,
# we're stuck using VTK5 for now, which has a cmake config
# that's hard to work with.

if(NOT VTK_DIR)
  set(VTK_DIR ${FS_PACKAGES_DIR}/vtk/5.6)
endif()

find_package(VTK HINTS ${VTK_DIR} NO_MODULE)

if(VTK_FOUND)

  # overwrite VTK_LIBRARIES with the absolute paths
  library_paths(
    NAME VTK_LIBRARIES
    LIBDIR ${VTK_LIBRARY_DIRS}
    LIBRARIES
    vtkverdict
    vtkGraphics
    vtkexpat
    vtkfreetype
    vtktiff
    vtkjpeg
    vtkmetaio
    vtkpng
    vtkzlib
    vtksqlite
    vtkImaging
    vtkFiltering
    vtkCommon
    vtksys
    vtkGenericFiltering
    vtkexoIIc
    vtkNetCDF
    vtkVolumeRendering
    vtkRendering
    vtkftgl
    vtkWidgets
    vtkHybrid
    vtkIO
    vtkDICOMParser
  )

  library_paths(
    NAME VTK_TCL_LIBRARIES
    LIBDIR ${VTK_LIBRARY_DIRS}
    LIBRARIES
    vtkImagingTCL
    vtkVolumeRenderingTCL
    vtkRenderingTCL
    vtkFilteringTCL
    vtkWidgetsTCL
    vtkHybridTCL
    vtkGraphicsTCL
    vtkImagingTCL
    vtkIOTCL
    vtkCommonTCL
  )

  # vtkWrapTcl and vtkWrapTclInit commands (required for qdec and vtkutils)
  if(NOT VTK_WRAP_TCL_EXE OR NOT VTK_WRAP_TCL_INIT_EXE)
    message(FATAL_ERROR "VTK must be built with VTK_WRAP_TCL ON")
  endif()

  # create a simple cmake function to use vtkWrapTcl
  function(vtk_wrap_tcl INFILE OUTFILE)
    add_custom_command(
      OUTPUT  ${OUTFILE}
      DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/${INFILE}
      COMMAND ${VTK_WRAP_TCL_EXE} ${CMAKE_CURRENT_SOURCE_DIR}/${INFILE}
              ${VTK_LIBRARY_DIRS}/hints 1 ${OUTFILE})
  endfunction()

  # create a simple cmake function to use vtkWrapTclInit
  function(vtk_wrap_tcl_init INFILE OUTFILE)
    add_custom_command(
      OUTPUT  ${OUTFILE}
      DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/${INFILE}
      COMMAND ${VTK_WRAP_TCL_INIT_EXE} ${CMAKE_CURRENT_SOURCE_DIR}/${INFILE} ${OUTFILE})
  endfunction()

  # make sure the vtk path gets linked from the lib directory during install
  if(VTK_INSTALL_PREFIX)
    symlink(${VTK_INSTALL_PREFIX} ${CMAKE_INSTALL_PREFIX}/lib/vtk)
  endif()

endif()
