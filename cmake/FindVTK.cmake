# VTK Find Module

# VTK6+ has a useful cmake config that sets up
# variables and component paths nicely. Unfortunately,
# we're stuck using VTK5 for now, which has a cmake config
# that's hard to work with.

if(NOT VTK_DIR)
  if(EXISTS ${FS_PACKAGES_DIR}/vtk/5.10.1)
    set(VTK_DIR ${FS_PACKAGES_DIR}/vtk/5.10.1)
  else()
    set(VTK_DIR ${FS_PACKAGES_DIR}/vtk/5.6.0)
  endif()
endif()

find_package(VTK HINTS ${VTK_DIR} NO_MODULE)

if(VTK_FOUND)

  # these are the libraries we need
  set(LIBS
    vtkverdict
    vtkGraphics
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
  )

  foreach(LIB ${LIBS})
    set(VTK_LIBRARIES ${VTK_LIBRARIES} ${VTK_LIBRARY_DIRS}/lib${LIB}.a)
  endforeach()

  # vtkWrapTcl command (required for the vtkutils)
  if(NOT EXISTS ${VTK_WRAP_TCL_EXE})
    message(FATAL_ERROR "VTK must be built with VTK_WRAP_TCL ON")
  endif()

  # make sure the vtk path gets linked from the lib directory during install
  symlink(${VTK_INSTALL_PREFIX} ${CMAKE_INSTALL_PREFIX}/lib/vtk)

endif()
