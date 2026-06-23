# FLTK is required
find_package(FLTK REQUIRED)
include_directories(${FLTK_INCLUDE_DIR})

# VTK as well
if(NOT VTK_FOUND)
  if(HOST_OS MATCHES "Rocky9" OR HOST_OS_MATCHES "Ubuntu24")
     find_package(VTK REQUIRED COMPONENTS CommonCore)
     include_directories(${VTK_INCLUDE_DIRS})
  else()
     find_package(VTK REQUIRED)
     include(${VTK_USE_FILE})
  endif()
endif()

# FLTK's fluid generates C++ files in the build directory
include_directories(${GEMS_SOURCE_DIR}/GUI)
include_directories(${GEMS_BINARY_DIR}/GUI)

# add libary
add_library(kvlGEMSGUI vtkFlRenderWindowInteractor.cxx kvlImageViewer.cxx)
target_link_libraries(kvlGEMSGUI
  ${ITK_LIBRARIES}
  ${VTK_LIBRARIES}
  ${FLTK_LIBRARIES}
)

if(LINK_SHARED_AS_NEEDED)
   target_link_options(kvlGEMSGUI PRIVATE "LINKER:--as-needed")
endif()

# add an executable
fltk_wrap_ui(kvlEstimateAtlasParameters kvlAtlasParameterEstimationConsoleGUI.fl)

add_executable(kvlEstimateAtlasParameters
  ${kvlEstimateAtlasParameters_FLTK_UI_SRCS}
  kvlEstimateAtlasParameters.cxx
  kvlAtlasParameterEstimationConsole.cxx
)
target_link_libraries(kvlEstimateAtlasParameters kvlGEMSCommon_dynmesh kvlGEMSGUI)

# add executable
fltk_wrap_ui(kvlViewMeshCollectionWithGUI kvlAtlasMeshViewingConsoleGUI.fl)

add_executable(kvlViewMeshCollectionWithGUI
  ${kvlViewMeshCollectionWithGUI_FLTK_UI_SRCS}
  kvlViewMeshCollectionWithGUI.cxx
  kvlAtlasMeshViewingConsole.cxx
)
target_link_libraries(kvlViewMeshCollectionWithGUI kvlGEMSCommon_dynmesh kvlGEMSGUI)
