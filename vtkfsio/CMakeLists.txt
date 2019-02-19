project(vtkfsio)

if(VTK_FOUND)

  include_directories(${FS_INCLUDE_DIRS} ${VTK_INCLUDE_DIRS})

  set(SOURCES
    vtkFSSurfaceLabelSource.cxx
    vtkFSSurfaceWriter.cxx
    vtkFSSurfaceScalarsReader.cxx
    vtkFSSurfaceSource.cxx
    vtkFSVolumeSource.cxx
  )

  if(SUPPRESS_WARNINGS)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-deprecated")
  endif()

  add_library(vtkfsio STATIC ${SOURCES})

endif()
