project(mris_pmake)

include_directories(${FS_INCLUDE_DIRS} ${CMAKE_CURRENT_SOURCE_DIR})

if(OpenCL_FOUND)
  include_directories(${OpenCL_INCLUDE_DIR})
endif()

set(SOURCES
  abs2rel.cpp
  asynch.cpp
  c_label.cpp
  C_mpmOverlay.cpp
  C_mpmProg.cpp
  c_SMessage.cpp
  c_SSocket.cpp
  c_surface.cpp
  c_vertex.cpp
  dijkstra.cpp
  env.cpp
  general.cpp
  help.cpp
  legacy.cpp
  mris_pmake.cpp
  oclDijkstraKernel.cpp
  rel2abs.cpp
  scanopt.cpp
  surface.cpp
)

if(OpenCL_FOUND)
  set(SOURCES ${SOURCES} oclCommon.cpp)
endif()

add_executable(mris_pmake ${SOURCES})
add_help(mris_pmake mris_pmake.help.xml)
target_link_libraries(mris_pmake utils)

if(OpenCL_FOUND)
  target_link_libraries(mris_pmake ${OpenCL_LIBRARY})
endif()

install(TARGETS mris_pmake DESTINATION bin)
