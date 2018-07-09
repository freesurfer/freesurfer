
# if the VTK_DIR path is not supplied and we're at Martinos,
# use this default path on pubsw
set(MARTINOS_VTK_DIR /usr/pubsw/packages/vtk/5.6.0_sshpatch)
if(NOT VTK_DIR AND EXISTS ${MARTINOS_VTK_DIR})
  set(VTK_DIR ${MARTINOS_VTK_DIR})
endif()

# if VTK_DIR was not supplied and we're building otuside of Martinos,
# download and compile VTK as an external project
if(NOT VTK_DIR)
  include(ExternalProject)
  # todo external package...
endif()

# now, actually run the find_package to set VTK variables
find_package(VTK PATHS ${VTK_DIR} NO_DEFAULT_PATH NO_MODULE REQUIRED)

# get the libraries
foreach(lib vtkFiltering vtkHybrid vtkRendering vtkGraphics vtkImaging vtkftgl vtkIO vtkCommon vtkDICOMParser vtksys vtkexoIIc vtkNetCDF)
  find_library(${lib}_LIBRARY PATHS ${VTK_LIBRARY_DIRS} NAMES ${lib})
  set(VTK_LIBRARIES ${VTK_LIBRARIES} ${${lib}_LIBRARY})
endforeach()
find_package_handle_standard_args(VTK DEFAULT_MSG VTK_LIBRARIES)
