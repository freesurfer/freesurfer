project(mris_pmake)

include_directories(${CMAKE_SOURCE_DIR}/include)

set(SOURCES
  abs2rel.c
  asynch.cpp
  c_label.cpp
  C_mpmProg.cpp
  c_SMessage.cpp
  c_SSocket.cpp
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
target_link_libraries(mris_pmake utils)

install(TARGETS mris_pmake DESTINATION bin)
