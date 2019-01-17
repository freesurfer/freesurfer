project(tetgen)

include_directories(${CMAKE_CURRENT_SOURCE_DIR})

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-array-bounds")

add_definitions(-DTETLIBRARY)

add_library(tetgen STATIC tetgen.cxx predicates.cxx)
