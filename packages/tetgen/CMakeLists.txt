project(tetgen)

add_definitions(-DTETLIBRARY)

add_library(tetgen STATIC tetgen.cxx predicates.cxx)
