project(hipsstubs)

include_directories(${FS_INCLUDE_DIRS})

SET(SOURCES hipsrepl.c hipsstubs.c)

add_library(hipsstubs STATIC ${SOURCES})
