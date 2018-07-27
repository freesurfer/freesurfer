project(log)

include_directories(${FS_INCLUDE_DIRS})

SET(SOURCES
  congraph.c
  connect.c
  h_logz.c
  lmedian.c
  map.c
  runfuncs.c
)

add_library(log STATIC ${SOURCES})
