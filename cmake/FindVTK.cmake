# VTK Find Module

# VTK6+ has a useful cmake config that sets up
# variables and component paths nicely. Unfortunately,
# we're stuck using VTK5 for now, which has a cmake config
# that's hard to work with.

if(NOT VTK_DIR)
  if(EXISTS ${FS_PACKAGES_DIR}/vtk/5.10.1)
    set(VTK_DIR ${FS_PACKAGES_DIR}/vtk/5.10.1)
  else()
    set(VTK_DIR ${FS_PACKAGES_DIR}/vtk/5.6)
  endif()
endif()

find_package(VTK HINTS ${VTK_DIR} NO_MODULE)

if(VTK_FOUND)

  # The order of the libraries is important
    
  # overwrite VTK_LIBRARIES with the absolute paths
  library_paths(
    NAME VTK_LIBRARIES
    LIBDIR ${VTK_LIBRARY_DIRS}
    LIBRARIES
    #
    vtkverdict
    vtkmetaio
    vtksqlite
    vtkexoIIc
    vtkNetCDF
    vtkNetCDF_cxx
    vtkDICOMParser
    vtkhdf5
    vtkhdf5_hl
    LSDyna
    #
    vtkWidgets
    vtkHybrid
    vtkVolumeRendering
    vtkRendering
    vtkIO
    vtkGenericFiltering
    vtkGraphics
    vtkImaging
    vtkFiltering
    #
    vtkftgl
    vtktiff
    vtkjpeg
    vtkpng
    vtkzlib
    vtkexpat
    vtkfreetype
    #
    vtkCommon
    vtksys
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

  if(NOT APPLE)
    # install the shared libraries to the freesurfer lib directory
    file(GLOB VTK_LIBRARIES_TO_INSTALL "${VTK_LIBRARY_DIRS}/lib*.so*")
    if(VTK_LIBRARIES_TO_INSTALL)
      install(PROGRAMS ${VTK_LIBRARIES_TO_INSTALL} DESTINATION lib/vtk)
      # add vtk library directory to rpath
      set(CMAKE_INSTALL_RPATH "$ORIGIN/../lib/vtk:${CMAKE_INSTALL_RPATH}")
    endif()
  endif()

endif()
