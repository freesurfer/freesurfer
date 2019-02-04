project(mris_pmake)

include_directories(
  ${FS_INCLUDE_DIRS}
  ${CMAKE_SOURCE_DIR}/opencl_algorithms
  ${CMAKE_CURRENT_SOURCE_DIR}
)

if(OpenCL_FOUND)
  include_directories(${OpenCL_INCLUDE_DIR})
endif()

set(SOURCES
  abs2rel.c
  asynch.cpp
  c_SMessage.cpp
  c_SSocket.cpp
  c_label.cpp
  C_mpmProg.cpp
  c_surface.cpp
  c_vertex.cpp
  dijkstra.cpp
  env.cpp
  general.cpp
  help.cpp
  mris_pmake.cpp
  rel2abs.c
  scanopt.cpp
  surface.cpp
  C_mpmOverlay.cpp
  legacy.cpp
)


add_executable(mris_pmake ${SOURCES})
add_help(mris_pmake mris_pmake.help.xml)
target_link_libraries(mris_pmake utils opencl_algorithms)

install(TARGETS mris_pmake DESTINATION bin)
