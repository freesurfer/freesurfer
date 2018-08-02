project(hipsstubs)

include_directories(${FS_INCLUDE_DIRS})

set(SOURCES hipsrepl.c hipsstubs.c)

add_library(hipsstubs STATIC ${SOURCES})
