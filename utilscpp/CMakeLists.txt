project(utilscpp)

include_directories(${FS_INCLUDE_DIRS})

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++0x")

set(SOURCES
  argparse.cpp
  log.cpp
  mris_fastmarching.cpp
  globals.cpp
  face.cpp
  vertex.cpp
  loop.cpp
  fastloop.cpp
  surface.cpp
  patchdisk.cpp
  mris_topology.cpp
  segment.cpp
  patchdisk.cpp
)

add_library(utilscpp STATIC ${SOURCES})
target_link_libraries(utilscpp utils)
