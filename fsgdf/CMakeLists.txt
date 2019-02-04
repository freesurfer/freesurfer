project(fsgdf)

include_directories(${FS_INCLUDE_DIRS} ${X11_INCLUDE_DIR})

add_library(fsgdf STATIC fsgdf.c)
target_link_libraries(fsgdf utils)
