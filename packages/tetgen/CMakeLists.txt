project(tetgen)

include_directories(${CMAKE_CURRENT_SOURCE_DIR})

add_definitions(-DTETLIBRARY)

add_library(tetgen STATIC tetgen.cxx predicates.cxx)
