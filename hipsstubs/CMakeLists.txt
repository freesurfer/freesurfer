project(hipsstubs)

include_directories(${CMAKE_SOURCE_DIR}/include)

SET(SOURCES hipsrepl.c hipsstubs.c)

add_library(hipsstubs STATIC ${SOURCES})
