project(mris_topo_fixer)

include_directories(${FS_INCLUDE_DIRS})

set(SOURCES
  face.cpp
  fastloop.cpp
  globals.cpp
  loop.cpp
  mris_topo_fixer.cpp
  mris_topology.cpp
  patchdisk.cpp
  segment.cpp
  surface.cpp
  vertex.cpp
)

add_executable(mris_topo_fixer ${SOURCES})
target_link_libraries(mris_topo_fixer utils)
install(TARGETS mris_topo_fixer DESTINATION bin)
