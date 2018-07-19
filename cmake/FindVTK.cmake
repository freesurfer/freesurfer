# VTK Find Module

# VTK6+ has a useful cmake config that sets up
# variables and component paths nicely. Unfortunately,
# we're stuck using VTK5 for now, which has a cmake config
# that's hard to work with. So instead, we'll use our own find module

if(NOT VTK_DIR)
  set(VTK_DIR ${FS_PACKAGES_DIR}/vtk/5.6.0_sshpatch)
endif()

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

find_path(VTK_INCLUDE_DIR NAMES vtkIOStream.h PATHS ${VTK_DIR} PATH_SUFFIXES include/vtk-5.6)
find_path(VTK_LIBRARY_DIR NAMES libvtkIO.so PATHS ${VTK_DIR} PATH_SUFFIXES lib/vtk-5.6)

# search for the libraries
foreach(LIB ${LIBS})
  find_library(tmp PATHS ${VTK_LIBRARY_DIR} NAMES ${LIB})
  set(VTK_LIBRARIES ${VTK_LIBRARIES} ${tmp})
  unset(tmp CACHE)  # this is necessary for find_library to work (plus it clears it from the cache)
endforeach()

find_package_handle_standard_args(VTK DEFAULT_MSG VTK_INCLUDE_DIR VTK_LIBRARIES)

# vtkWrapTcl command (required for the vtkutils)
if(EXISTS ${VTK_DIR}/bin/vtkWrapTcl)
  set(VTK_WRAP_TCL ${VTK_DIR}/bin/vtkWrapTcl)
endif()

# make sure the vtk path gets linked from the lib directory during install
symlink(${VTK_DIR} ${CMAKE_INSTALL_PREFIX}/lib/vtk)
